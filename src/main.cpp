/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file main.cpp
 * @brief Main entry point for the 3D incompressible CFD solver
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <iostream>
#include <iomanip>
#include <chrono>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"
#include "CFDApplication.h"
#include "Comm.h"
#include "ErrorHandler.h"
#include "Logger.h"
#include "PETScRuntime.h"

// *********************************** main ***********************************

int main(int argc, char* argv[])
{
    // Start timing the total execution
    const auto startTime = std::chrono::high_resolution_clock::now();

    // PETSc/MPI runtime initialization
    const PETScRuntime petscRuntime;

    // One console, whole-run failure: non-master ranks print nothing, and
    // a fatal error on any rank aborts every rank. Single-rank runs keep
    // std::abort — same exit code and core dump as a serial program
    Logger::init(Comm::master());

    if (Comm::parallelRun())
    {
        setAbortHandler(&Comm::abortAllRanks);
    }

    // Scaffolding until the mesh is decomposed across ranks (Phase 4/5)
    if (Comm::parallelRun())
    {
        FatalError
        (
            "Multi-rank runs are not supported yet: the mesh is not "
            "decomposed. Run on a single rank."
        );
    }

    std::cout << R"(
  ~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~··~·~

  ████████╗██╗   ██╗██████╗ ██████╗ ██╗  ██╗   ██╗███████╗███████╗
  ╚══██╔══╝██║   ██║██╔══██╗██╔══██╗██║  ╚██╗ ██╔╝╚══███╔╝██╔════╝
     ██║   ██║   ██║██████╔╝██████╔╝██║   ╚████╔╝   ███╔╝ █████╗
     ██║   ██║   ██║██╔══██╗██╔══██╗██║    ╚██╔╝   ███╔╝  ██╔══╝
     ██║   ╚██████╔╝██║  ██║██████╔╝███████╗██║   ███████╗███████╗
     ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚══════╝╚═╝   ╚══════╝╚══════╝

           3D Incompressible Navier-Stokes Solver v1.0

  ~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~·~··~·~
)" << '\n';

    std::cout
        << "Running with precision: " << SCALAR_MODE << '\n';

    std::cout 
        << std::fixed << std::setprecision(6);

    FilePath caseFile = "../defaultCase";

    if (argc > 1)
    {
        caseFile = argv[1];

        std::cout
            << "Using case file: " << caseFile << '\n';
    }
    else
    {
        std::cout
            << "Using default case: " << caseFile << '\n';
    }

    CFDApplication app(caseFile);
    app.run();

    const auto endTime = std::chrono::high_resolution_clock::now();

    const auto duration =
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

    std::cout << '\n';
    Logger::sectionHeader("Simulation Complete");
    Logger::keyValue
    (
        "Total execution time",
        std::to_string(duration.count()) + " seconds"
    );
    Logger::iterationFooter();
    std::cout << '\n';

    return 0;
}
