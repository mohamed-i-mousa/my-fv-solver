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
#include <cmath>
#include <iostream>

// External library headers
#include <eigen3/Eigen/Core>
#include <omp.h>

// Project headers
#include "BoundaryConditionsLoader.h"
#include "BoundaryConditions.h"
#include "CaseConfiguration.h"
#include "CaseReader.h"
#include "Forces.h"
#include "Logger.h"
#include "MeshCreator.h"
#include "PostProcess.h"
#include "SolverSetup.h"
#include "MomentumTransport.h"
#include "PvdTimeSeries.h"

// ***************************** Internal Helpers *****************************

namespace
{

void initParallelism(int numThreads)
{
    omp_set_num_threads(numThreads);
    Eigen::setNbThreads(numThreads);

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

    // Set up the PVD time series and the optional force-history CSV
    const FilePath pvdFile = PostProcess::pvdPathFor(config);
    VTK::writePVDTimeSeriesHeader(pvdFile);

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

    // Export the initial condition (t = 0) before the transient solve
    PostProcess::exportTimeStep
    (
        pvdFile,
        S(0.0),
        0,
        mesh,
        *modules.solver,
        *modules.turbulenceModel,
        config
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
            PostProcess::exportTimeStep
            (
                pvdFile,
                time,
                step,
                mesh,
                *modules.solver,
                *modules.turbulenceModel,
                config
            );
        }
    }

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

    // Create mesh
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
