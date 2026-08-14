/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LeastSquares.cpp
 * @brief Implementation of the weighted least-squares gradient scheme
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "LeastSquares.h"

// External library headers
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>

// Project headers
#include "ErrorHandler.h"

// ************************* Special Member Functions *************************

LeastSquares::LeastSquares
(
    const Mesh& mesh,
    const BoundaryConditions& bc
)
:
    GradientScheme(mesh, bc)
{
    precomputeInverseATA();
}

// ****************************** Public Methods ******************************

Vector LeastSquares::cellGradient
(
    Field field,
    const ScalarField& phi,
    Index cellIndex
) const
{
    const Cell& cell = mesh().cells()[cellIndex];

    Scalar b0 = S(0.0);
    Scalar b1 = S(0.0);
    Scalar b2 = S(0.0);

    // Part 1: Internal neighbor cells contribution to ATb
    for (Index neighborIdx : cell.neighborCellIndices())
    {
        const Cell& neighbor = mesh().cells()[neighborIdx];
        const Vector r = neighbor.centroid() - cell.centroid();

        const Scalar rMagSqr = magnitudeSquared(r);
        const Scalar w = S(1.0) / (rMagSqr + smallValue);

        const Scalar wDeltaPhi =
            w * (phi[neighborIdx] - phi[cellIndex]);

        b0 += wDeltaPhi * r.x();
        b1 += wDeltaPhi * r.y();
        b2 += wDeltaPhi * r.z();
    }

    // Part 2: Boundary faces contribution to ATb
    for (Index faceIdx : cell.faceIndices())
    {
        const Face& face = mesh().faces()[faceIdx];

        if (!face.isBoundary()) continue;

        const Vector r = face.centroid() - cell.centroid();
        const Scalar rMagSqr = magnitudeSquared(r);
        const Scalar w = S(1.0) / (rMagSqr + smallValue);

        const Index bIdx = bcManager().boundaryIdx(face.idx());
        const Scalar phiBoundary =
            bcManager().boundaryType(field, bIdx).faceValue
            (
                phi[face.ownerCell()],
                bcManager().normalDistance(bIdx),
                bcManager().normal(bIdx),
                bcManager().ownerVelocity(bIdx)
            );

        const Scalar wDeltaPhi =
            w * (phiBoundary - phi[cellIndex]);

        b0 += wDeltaPhi * r.x();
        b1 += wDeltaPhi * r.y();
        b2 += wDeltaPhi * r.z();
    }

    // g = inv(ATA) * ATb  (symmetric 3x3 mat-vec multiply)
    const auto& inv = invATA_[cellIndex];

    return Vector
    (
        inv[0]*b0 + inv[1]*b1 + inv[2]*b2,
        inv[1]*b0 + inv[3]*b1 + inv[4]*b2,
        inv[2]*b0 + inv[4]*b1 + inv[5]*b2
    );
}


// ****************************** Private Methods *****************************

void LeastSquares::precomputeInverseATA()
{
    // Aliases for Eigen types
    using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using CholeskySolver = Eigen::LLT<Matrix3>;
    using LUSolver = Eigen::FullPivLU<Matrix3>;

    // Sized over every cell; ghosts are not gradient sites and stay zero
    invATA_.resize(mesh().numCells());

    const Count numOwnedCells = mesh().numOwnedCells();

    Count degenerateCells = 0;

    for (Index cellIdx = 0; cellIdx < numOwnedCells; ++cellIdx)
    {
        Matrix3 ATA;
        Vector3 rVector;
        const Cell& cell = mesh().cells()[cellIdx];

        ATA.setZero();

        // Neighbor cells contribution (purely geometric)
        for (Index neighborIdx : cell.neighborCellIndices())
        {
            const Vector r =
                mesh().cells()[neighborIdx].centroid()
              - cell.centroid();

            const Scalar rMagSqr = magnitudeSquared(r);
            const Scalar w = S(1.0) / (rMagSqr + smallValue);

            rVector << r.x(), r.y(), r.z();
            ATA.noalias() += w * (rVector * rVector.transpose());
        }

        // Boundary faces contribution (purely geometric)
        for (Index faceIdx : cell.faceIndices())
        {
            const Face& face = mesh().faces()[faceIdx];

            if (!face.isBoundary()) continue;

            const Vector r = face.centroid() - cell.centroid();
            const Scalar rMagSqr = magnitudeSquared(r);
            const Scalar w = S(1.0) / (rMagSqr + smallValue);

            rVector << r.x(), r.y(), r.z();
            ATA.noalias() += w * (rVector * rVector.transpose());
        }

        // Invert ATA and store symmetric result
        Matrix3 inv;
        bool inverted = false;

        CholeskySolver llt(ATA);

        if (llt.info() == Eigen::Success)
        {
            inv = llt.solve(Matrix3::Identity());
            inverted = true;
        }
        else
        {
            LUSolver lu(ATA);

            if (lu.isInvertible())
            {
                inv = lu.inverse();
                inverted = true;
            }
        }

        if (inverted)
        {
            // Store upper triangle: {xx, xy, xz, yy, yz, zz}
            invATA_[cellIdx] =
            {
                inv(0,0), inv(0,1), inv(0,2),
                          inv(1,1), inv(1,2),
                                    inv(2,2)
            };
        }
        else
        {
            invATA_[cellIdx] = {0, 0, 0, 0, 0, 0};
            ++degenerateCells;
        }
    }

    if (degenerateCells > 0)
    {
        Warning
        (
            std::to_string(degenerateCells)
          + " cells have degenerate least-squares"
            " matrices (gradient will be zero)"
        );
    }
}
