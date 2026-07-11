/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Forces.cpp
 * @brief Aerodynamic force calculation on a wall patch
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Forces.h"

// Standard library headers
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

// Project headers
#include "Vector.h"
#include "BoundaryConditions.h"
#include "BoundaryPatch.h"
#include "CaseConfiguration.h"
#include "ErrorHandler.h"
#include "Face.h"
#include "Field.h"
#include "Logger.h"
#include "Mesh.h"
#include "MomentumTransport.h"
#include "Comm.h"
#include "Reduce.h"
#include "TurbulenceModel.h"

// ***************************** namespace Forces ****************************

namespace Forces
{

// ***************************** Internal Helpers *****************************

namespace
{

// Drag/lift loads split into pressure and friction contributions
struct AeroForces
{
    Scalar pressureDrag;
    Scalar frictionDrag;
    Scalar pressureLift;
    Scalar frictionLift;
};

// Derive a forces output-file path from the configured output file name
FilePath forcesFilePath
(
    const FilePath& outputFilename,
    const FilePath& suffix
)
{
    static const std::array<FilePath, 4> extensions =
    {
        FilePath{".vtkhdf"},
        FilePath{".vtu"},
        FilePath{".vtp"},
        FilePath{".pvd"}
    };

    FilePath path = outputFilename;

    for (const FilePath& extension : extensions)
    {
        if (path.ends_with(extension))
        {
            path.resize(path.size() - extension.size());
            break;
        }
    }

    return path + suffix;
}

// Dynamic load 0.5 * rho * Vref^2 * A used to non-dimensionalize forces
Scalar referenceDynamicLoad(const CaseConfiguration& config)
{
    const Scalar referenceVelocity = magnitude(config.referenceVelocity);
    const Scalar dynamicLoad =
        S(0.5) * config.rho * referenceVelocity * referenceVelocity
      * config.referenceArea;

    if (dynamicLoad <= vSmallValue)
    {
        FatalError
        (
            "Reference dynamic load is zero; force coefficients are "
            "undefined. Check that density, referenceVelocity, and "
            "referenceArea are all non-zero."
        );
    }

    return dynamicLoad;
}

AeroForces computeForces
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    // The forces patch may be absent from this rank's submesh: found on
    // no rank at all is a case error, found elsewhere means this rank
    // integrates nothing and only joins the collective sum below
    const BoundaryPatch* patch = nullptr;

    for (const BoundaryPatch& meshPatch : mesh.patches())
    {
        if (meshPatch.patchName() == config.forcesPatch)
        {
            patch = &meshPatch;
            break;
        }
    }

    if (globalSum(Count{patch != nullptr ? 1u : 0u}) == 0)
    {
        FatalError
        (
            "Forces patch '" + config.forcesPatch
          + "' not found on any rank"
        );
    }

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();
    const FaceData<Scalar> wallShearStress =
        turbulence.wallShearStress(Ux, Uy, Uz);

    const FaceListRef faces = mesh.faces();

    Scalar pressureForceX = S(0.0);
    Scalar pressureForceY = S(0.0);
    Scalar pressureForceZ = S(0.0);
    Scalar frictionForceX = S(0.0);
    Scalar frictionForceY = S(0.0);
    Scalar frictionForceZ = S(0.0);

    const Index firstFaceIdx = patch ? patch->firstFaceIdx() : 1;
    const Index lastFaceIdx = patch ? patch->lastFaceIdx() : 0;

    #pragma omp parallel for schedule(static) \
        reduction(+:pressureForceX, pressureForceY, pressureForceZ) \
        reduction(+:frictionForceX, frictionForceY, frictionForceZ)
    for
    (
        Index faceIdx = firstFaceIdx;
        faceIdx <= lastFaceIdx;
        ++faceIdx
    )
    {
        const Face& face = faces[faceIdx];
        const Index cellIdx = face.ownerCell();
        const Vector& normal = face.normal();

        // Pressure force from the kinematic pressure
        const Scalar pressureFace =
            bcManager.boundaryFaceValue(Field::p, pressure, face);
        const Vector pressureContribution =
            (config.rho * pressureFace * face.projectedArea()) * normal;
        pressureForceX += pressureContribution.x();
        pressureForceY += pressureContribution.y();
        pressureForceZ += pressureContribution.z();

        // Skin-friction force
        const Vector cellVelocity(Ux[cellIdx], Uy[cellIdx], Uz[cellIdx]);
        const Scalar normalVelocity = dot(cellVelocity, normal);
        const Vector tangentVelocity =
            cellVelocity - normalVelocity * normal;
        const Scalar tangentMagnitude = magnitude(tangentVelocity);

        if (tangentMagnitude > vSmallValue)
        {
            const Vector shearDirection = tangentVelocity / tangentMagnitude;
            const Scalar shearStress =
                config.rho * wallShearStress[face.idx()];

            const Vector frictionContribution =
                (shearStress * face.contactArea()) * shearDirection;
            frictionForceX += frictionContribution.x();
            frictionForceY += frictionContribution.y();
            frictionForceZ += frictionContribution.z();
        }
    }

    // Each rank integrates its share of the patch; combine onto every rank
    const Vector pressureForce =
        globalSum(Vector(pressureForceX, pressureForceY, pressureForceZ));
    const Vector frictionForce =
        globalSum(Vector(frictionForceX, frictionForceY, frictionForceZ));

    const Vector& dragDir = config.dragDirection;
    const Vector& liftDir = config.liftDirection;

    return AeroForces
    {
        .pressureDrag = dot(pressureForce, dragDir),
        .frictionDrag = dot(frictionForce, dragDir),
        .pressureLift = dot(pressureForce, liftDir),
        .frictionLift = dot(frictionForce, liftDir)
    };
}

} // namespace

// ***************************** Force Reporting ******************************

void reportForces
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    const AeroForces forces =
        computeForces(solver, turbulence, mesh, bcManager, config);

    const Vector& dragDir = config.dragDirection;
    const Vector& liftDir = config.liftDirection;

    const Scalar pressureDrag = forces.pressureDrag;
    const Scalar frictionDrag = forces.frictionDrag;
    const Scalar totalDrag = pressureDrag + frictionDrag;

    const Scalar pressureLift = forces.pressureLift;
    const Scalar frictionLift = forces.frictionLift;
    const Scalar totalLift = pressureLift + frictionLift;

    // Non-dimensionalise: C = F / (0.5 * rho * Vref^2 * A).
    const Scalar dynamicLoad = referenceDynamicLoad(config);

    const Scalar pressureCd = pressureDrag / dynamicLoad;
    const Scalar frictionCd = frictionDrag / dynamicLoad;
    const Scalar totalCd = totalDrag / dynamicLoad;

    const Scalar pressureCl = pressureLift / dynamicLoad;
    const Scalar frictionCl = frictionLift / dynamicLoad;
    const Scalar totalCl = totalLift / dynamicLoad;

    // Write the breakdown to a text file beside the VTK output
    const FilePath outputPath =
        forcesFilePath(config.vtkOutputFilename, "_forces.txt");

    // Every rank reduced the forces above; only the master writes
    if (Comm::master())
    {
        std::ofstream file(outputPath);
        if (!file.is_open())
        {
            FatalError("Failed to open forces output file: " + outputPath);
        }

        file << std::scientific << std::setprecision(6);
        file
            << "Aerodynamic forces" << '\n'
            << "Patch          : " << config.forcesPatch << '\n'
            << "Drag direction : " << dragDir << '\n'
            << "Lift direction : " << liftDir << '\n'
            << "Reference U    : " << config.referenceVelocity << '\n'
            << '\n'
            << std::left << std::setw(12) << "Force"
            << std::right
            << std::setw(16) << "Pressure"
            << std::setw(16) << "Friction"
            << std::setw(16) << "Total" << '\n'
            << std::left << std::setw(12) << "Drag [N]"
            << std::right
            << std::setw(16) << pressureDrag
            << std::setw(16) << frictionDrag
            << std::setw(16) << totalDrag << '\n'
            << std::left << std::setw(12) << "Lift [N]"
            << std::right
            << std::setw(16) << pressureLift
            << std::setw(16) << frictionLift
            << std::setw(16) << totalLift << '\n'
            << '\n'
            << "Force coefficients (dimensionless), referenceArea = "
            << config.referenceArea << '\n'
            << std::left << std::setw(12) << "Coeff"
            << std::right
            << std::setw(16) << "Pressure"
            << std::setw(16) << "Friction"
            << std::setw(16) << "Total" << '\n'
            << std::left << std::setw(12) << "Cd"
            << std::right
            << std::setw(16) << pressureCd
            << std::setw(16) << frictionCd
            << std::setw(16) << totalCd << '\n'
            << std::left << std::setw(12) << "Cl"
            << std::right
            << std::setw(16) << pressureCl
            << std::setw(16) << frictionCl
            << std::setw(16) << totalCl << '\n';
        file.close();
    }

    // Print the breakdown to the console as two compact tables
    std::cout << '\n';
    Logger::sectionHeader("Aerodynamic Forces");
    Logger::subsection("Patch: " + config.forcesPatch);

    Logger::breakdownHeader("Forces [N]");
    Logger::breakdownRow("Drag", pressureDrag, frictionDrag, totalDrag);
    Logger::breakdownRow("Lift", pressureLift, frictionLift, totalLift);

    Logger::breakdownHeader("Coefficients");
    Logger::breakdownRow("Cd", pressureCd, frictionCd, totalCd);
    Logger::breakdownRow("Cl", pressureCl, frictionCl, totalCl);

    Logger::subsection("Output file: " + outputPath);
    Logger::iterationFooter();
}

// ************************* Transient Force History *************************

void writeForceHistoryHeader(const CaseConfiguration& config)
{
    const FilePath csvPath =
        forcesFilePath(config.vtkOutputFilename, "_forces.csv");

    if (!Comm::master())
    {
        return;
    }

    std::ofstream file(csvPath);
    if (!file.is_open())
    {
        FatalError("Failed to open forces history file: " + csvPath);
    }

    file << "time,pressureDrag,frictionDrag,totalDrag,"
         << "pressureLift,frictionLift,totalLift,Cd,Cl" << '\n';
}


void appendForceHistory
(
    Scalar time,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const CaseConfiguration& config
)
{
    const AeroForces forces =
        computeForces(solver, turbulence, mesh, bcManager, config);

    const Scalar totalDrag = forces.pressureDrag + forces.frictionDrag;
    const Scalar totalLift = forces.pressureLift + forces.frictionLift;

    const Scalar dynamicLoad = referenceDynamicLoad(config);

    const Scalar totalCd = totalDrag / dynamicLoad;
    const Scalar totalCl = totalLift / dynamicLoad;

    // Every rank reduced the forces above; only the master writes
    if (!Comm::master())
    {
        return;
    }

    const FilePath csvPath =
        forcesFilePath(config.vtkOutputFilename, "_forces.csv");

    std::ofstream file(csvPath, std::ios::app);
    if (!file.is_open())
    {
        FatalError("Failed to open forces history file: " + csvPath);
    }

    file << std::scientific << std::setprecision(6)
         << time << ','
         << forces.pressureDrag << ',' << forces.frictionDrag << ','
         << totalDrag << ','
         << forces.pressureLift << ',' << forces.frictionLift << ','
         << totalLift << ','
         << totalCd << ',' << totalCl << '\n';
}

} // namespace Forces
