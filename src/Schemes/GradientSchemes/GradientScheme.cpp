/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GradientScheme.cpp
 * @brief Implementation of scheme-independent gradient services
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "GradientScheme.h"

// Standard library headers
#include <algorithm>

// Project headers
#include "ErrorHandler.h"
#include "LeastSquares.h"
#include "RuntimeSelection.h"

// ************************* Special Member Functions *************************

GradientScheme::GradientScheme
(
    const Mesh& mesh,
    const BoundaryConditions& bc
)
:
    mesh_{mesh},
    bcManager_{bc}
{}

// **************************** Runtime Selection ****************************

std::unique_ptr<GradientScheme> GradientScheme::create
(
    const Name& schemeName,
    const Mesh& mesh,
    const BoundaryConditions& bc
)
{
    if (schemeName == "leastSquares")
    {
        return std::make_unique<LeastSquares>(mesh, bc);
    }

    RuntimeSelection::unknownSelection
    (
        "gradient scheme",
        schemeName,
        availableSchemes()
    );
}


NameList GradientScheme::availableSchemes()
{
    return {"leastSquares"};
}

// ****************************** Public Methods ******************************

Vector GradientScheme::faceGradient
(
    Field field,
    const ScalarField& phi,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    Index faceIndex
) const
{
    const Face& f = mesh_.faces()[faceIndex];
    const Index P = f.ownerCell();

    if (f.isBoundary())
    {
        return
            boundaryFaceGradient
            (
                field,
                phi,
                gradPhiP,
                f
            );
    }

    const Index N = f.neighborCell().value();
    const Vector dPN =
        mesh_.cells()[N].centroid() - mesh_.cells()[P].centroid();
    const Scalar dPNMag = magnitude(dPN);
    const Vector ePN = dPN / (dPNMag + vSmallValue);
    const Vector gradAvg = averageFaceGradient(f, gradPhiP, gradPhiN);
    const Scalar phiDiff = phi[N] - phi[P];
    const Scalar correction =
        (phiDiff / (dPNMag + vSmallValue))
      - dot(gradAvg, ePN);

    return gradAvg + correction * ePN;
}


void GradientScheme::limitGradient
(
    Field field,
    const ScalarField& phi,
    VectorField& gradPhi
) const
{
    const Count numCells = mesh_.numOwnedCells();

    // Symmetry faces excluded: mirroring needs all three components
    const bool isVelocity =
        field == Field::Ux || field == Field::Uy || field == Field::Uz;

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Cell& cell = mesh_.cells()[cellIdx];

        // Find min/max among cell P and all its neighbors
        Scalar phiMin = phi[cellIdx];
        Scalar phiMax = phi[cellIdx];

        for (Index neighborIdx : cell.neighborCellIndices())
        {
            phiMin = std::min(phiMin, phi[neighborIdx]);
            phiMax = std::max(phiMax, phi[neighborIdx]);
        }

        // Boundary values included so the limiter does not clip k at walls
        for (Index faceIdx : cell.faceIndices())
        {
            const Face& f = mesh_.faces()[faceIdx];
            if (!f.isBoundary()) continue;

            if (isVelocity && bcManager_.excludedFromVelocityHull(f.idx()))
            {
                continue;
            }

            const Index bIdx = bcManager_.boundaryIdx(f.idx());
            const Scalar phiBound =
                bcManager_.boundaryType(field, bIdx).faceValue
                (
                    phi[f.ownerCell()],
                    bcManager_.normalDistance(bIdx),
                    bcManager_.normal(bIdx),
                    bcManager_.ownerVelocity(bIdx)
                );
            phiMin = std::min(phiMin, phiBound);
            phiMax = std::max(phiMax, phiBound);
        }

        // Compute Barth-Jespersen limiter: min alpha over all faces
        Scalar alpha = S(1.0);

        for (Index faceIdx : cell.faceIndices())
        {
            const Face& f = mesh_.faces()[faceIdx];

            if
            (
                isVelocity && f.isBoundary()
             && bcManager_.excludedFromVelocityHull(f.idx())
            )
            {
                continue;
            }

            const Vector r = f.centroid() - cell.centroid();
            const Scalar phiFace = phi[cellIdx] + dot(gradPhi[cellIdx], r);
            const Scalar delta = phiFace - phi[cellIdx];

            if (delta > smallValue)
            {
                alpha = std::min(alpha, (phiMax - phi[cellIdx]) / delta);
            }
            else if (delta < -smallValue)
            {
                alpha = std::min(alpha, (phiMin - phi[cellIdx]) / delta);
            }
        }

        alpha = std::clamp(alpha, S(0.0), S(1.0));

        gradPhi[cellIdx] = alpha * gradPhi[cellIdx];
    }
}


void GradientScheme::fieldGradient
(
    Field field,
    const ScalarField& phi,
    VectorField& gradPhi
) const
{
    const Count numCells = mesh_.numOwnedCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        gradPhi[cellIdx] = cellGradient(field, phi, cellIdx);
    }

    limitGradient(field, phi, gradPhi);
}

// ****************************** Private Methods *****************************

Vector GradientScheme::averageFaceGradient
(
    const Face& internalFace,
    const Vector& gradPhiP,
    const Vector& gradPhiN
) const
{
    const Scalar dPf = internalFace.dPfMag();
    const Scalar dNf = internalFace.dNfMag().value();
    const Scalar totalDist = dPf + dNf;

    const Scalar gP = dNf / (totalDist + vSmallValue);
    const Scalar gN = dPf / (totalDist + vSmallValue);

    return gP * gradPhiP + gN * gradPhiN;
}


Vector GradientScheme::boundaryFaceGradient
(
    Field field,
    const ScalarField& phi,
    const Vector& cellGradient,
    const Face& boundaryFace
) const
{
    const Vector tangentialGradient =
        cellGradient
      - dot(cellGradient, boundaryFace.normal())
      * boundaryFace.normal();

    // Normal gradient recovered from the boundary condition's face value
    const Index bIdx = bcManager_.boundaryIdx(boundaryFace.idx());
    const Scalar boundaryValue =
        bcManager_.boundaryType(field, bIdx).faceValue
        (
            phi[boundaryFace.ownerCell()],
            bcManager_.normalDistance(bIdx),
            bcManager_.normal(bIdx),
            bcManager_.ownerVelocity(bIdx)
        );
    const Scalar cellValue = phi[boundaryFace.ownerCell()];
    const Scalar dn = dot(boundaryFace.dPf(), boundaryFace.normal());
    const Scalar dPfMag = boundaryFace.dPfMag();

    // Stabilization: clamp dn to minNormalFraction_ * ||dPf||
    const Scalar dnStabilized =
        std::max(dn, minNormalFraction_ * dPfMag);

    const Scalar normalGradient =
        (boundaryValue - cellValue) / dnStabilized;

    return tangentialGradient + normalGradient * boundaryFace.normal();
}
