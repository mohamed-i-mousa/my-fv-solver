/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CFDApplication.cpp
 * @brief Top-level application orchestrator for the CFD solver
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "CFDApplication.h"

// Standard library headers
#include <optional>
#include <cmath>
#include <iostream>

// External library headers
#include <omp.h>

// Project headers
#include "BoundaryConditionsLoader.h"
#include "BoundaryConditions.h"
#include "CaseConfiguration.h"
#include "Comm.h"
#include "CaseReader.h"
#include "Forces.h"
#include "Logger.h"
#include "MeshCreator.h"
#include "PostProcess.h"
#include "SolverSetup.h"
#include "MomentumTransport.h"
#include "HDF5BoundaryData.h"
#include "HDF5CellData.h"

// ***************************** Internal Helpers *****************************

namespace
{

void initParallelism(int numThreads)
{
    omp_set_num_threads(numThreads);

    std::cout
        << "OpenMP threads: " << numThreads << '\n';
}


void runSteady
(
    SolverModules& modules,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    // Run the solver to convergence
    modules.solver->solve();

    // Post-process results
    PostProcess::reportStatistics(*modules.solver);

    // One shared VTKHDF file per grid; parallel runs write each rank's
    // owned cells as one piece of it, collectively via MPI-IO
    PostProcess::exportResults
    (
        *modules.solver,
        *modules.turbulenceModel,
        mesh,
        config
    );

    // Integrate aerodynamic forces on the configured wall patch
    if (config.forcesEnabled)
    {
        Forces::reportForces
        (
            *modules.solver,
            *modules.turbulenceModel,
            mesh,
            bcManager,
            config
        );
    }
}


void runTransient
(
    SolverModules& modules,
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    Logger::sectionHeader("Starting Transient Loop");

    if (config.forcesEnabled)
    {
        Forces::writeForceHistoryHeader(config);
    }

    Count numSteps = static_cast<Count>
    (
        std::floor(config.time.totalTime / config.time.timeStep + S(0.5))
    );
    if (numSteps == 0)
    {
        numSteps = 1;
    }

    // One shared VTKHDF file per grid (one piece per rank when
    // decomposed): geometry once, cell data appended per step
    VTK::HDF5CellData volumeWriter
    (
        PostProcess::cellDataPath(config),
        mesh,
        config.debug
    );
    VTK::HDF5BoundaryData boundaryWriter
    (
        PostProcess::boundaryDataPath(config),
        mesh,
        config.debug
    );

    volumeWriter.writeGeometry();
    boundaryWriter.writeGeometry();

    // Export the initial condition (t = 0) before the transient solve
    PostProcess::appendTimeStep
    (
        volumeWriter,
        boundaryWriter,
        S(0.0),
        *modules.solver,
        *modules.turbulenceModel
    );

    if (config.forcesEnabled)
    {
        Forces::appendForceHistory
        (
            S(0.0),
            mesh,
            bcManager,
            *modules.solver,
            *modules.turbulenceModel,
            config
        );
    }

    TransientFields prevStep;

    for (Count step = 1; step <= numSteps; ++step)
    {
        const Scalar time = S(step) * config.time.timeStep;

        modules.solver->solve(step, numSteps, time, &prevStep);

        if (config.forcesEnabled)
        {
            Forces::appendForceHistory
            (
                time,
                mesh,
                bcManager,
                *modules.solver,
                *modules.turbulenceModel,
                config
            );
        }

        const bool writeNow =
            (step % config.time.writingIntervals == 0) || (step == numSteps);

        if (writeNow)
        {
            PostProcess::appendTimeStep
            (
                volumeWriter,
                boundaryWriter,
                time,
                *modules.solver,
                *modules.turbulenceModel
            );
        }
    }

    volumeWriter.finalize();
    boundaryWriter.finalize();

    // Final reporting on the last time step's state
    PostProcess::reportStatistics(*modules.solver);

    if (config.forcesEnabled)
    {
        Forces::reportForces
        (
            *modules.solver,
            *modules.turbulenceModel,
            mesh,
            bcManager,
            config
        );
    }
}

} // namespace

// ************************* Special Member Functions *************************

CFDApplication::CFDApplication(const FilePath& caseFile)
:
    caseFile_{caseFile}
{}

CFDApplication::~CFDApplication() noexcept = default;

// ******************************** Solver Run ********************************

void CFDApplication::run()
{
    std::cout << '\n';
    Logger::sectionHeader("Loading Case");

    // Read the case file and load configuration
    CaseReader caseReader(caseFile_);
    const CaseConfiguration config = CaseConfig::loadConfiguration(caseReader);

    // Initialize parallelism
    initParallelism(static_cast<int>(config.numThreads));

    // Create mesh (decomposed across ranks in a parallel run)
    Mesh mesh = MeshCreator::create(config);

    // Load boundary conditions
    BoundaryConditions bcManager;
    BCLoader::load(caseReader, config, mesh, bcManager);

    // Configure solver
    SolverModules modules;
    SolverSetup::configure(modules, mesh, bcManager, config);
    SolverSetup::logSetup(modules, config);

    if (!modules.solver->isTransient())
    {
        runSteady(modules, mesh, bcManager, config);
    }
    else
    {
        runTransient(modules, mesh, bcManager, config);
    }
}
