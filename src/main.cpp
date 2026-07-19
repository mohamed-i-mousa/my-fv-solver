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
#include "MPIRuntime.h"
#include "PETScRuntime.h"

// *********************************** main ***********************************

int main(int argc, char* argv[])
{
    // Start timing the total execution
    const auto startTime = std::chrono::high_resolution_clock::now();

    // Initialize MPI and PETSc runtimes
    const MPIRuntime mpi;
    const PETScRuntime petsc;

    // Non-master processes will not output to std::cout
    if (!Comm::master()) { std::cout.setstate(std::ios::failbit); }

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
