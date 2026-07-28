/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TestMain.cpp
 * @brief Custom Catch2 main for the MPI + PETSc test executable
 *
 * @details Mirrors the runtime setup in src/main.cpp: an MPIRuntime (MPI_Init
 * plus Comm::init) then a PETScRuntime, both RAII.
 * A fatal error routes through FatalError's MPI_Abort, so no abort handler is
 * registered here. The per-rank failure count is summed so any rank's failure
 * fails the whole job.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <ios>
#include <iostream>
#include <string>

// External library headers
#include <catch2/catch_session.hpp>

// Project headers
#include "MPIRuntime.h"
#include "PETScRuntime.h"
#include "Comm.h"
#include "Reduce.h"
#include "Integer.h"
#include "StringTypes.h"

// ***************************** Internal Helpers *****************************

namespace
{

/// Reporter output file of a non-master rank
[[nodiscard]] FilePath rankReportFile(Index rank, Count size)
{
    return "mpiTestNp" + std::to_string(size)
         + "Rank" + std::to_string(rank) + ".log";
}

/// Construct the MPI and PETSc runtimes, run the session, and reduce the
/// per-rank failure count so any rank's failure fails the whole job
int runWithRuntimes(Catch::Session& session)
{
    const MPIRuntime mpi;
    const PETScRuntime petsc;

    const Index rank = Comm::myProcessorNum();
    const Count size = Comm::numProcessors();

    if (!Comm::master())
    {
        std::cout.setstate(std::ios::failbit);
        session.configData().defaultOutputFilename = rankReportFile(rank, size);
    }

    const int localFailures = session.run();

    if (!Comm::master() && localFailures > 0)
    {
        std::cout.clear();
        std::cout
            << "rank " << rank << ": " << localFailures
            << " failure(s), see " << rankReportFile(rank, size) << '\n';
    }

    const Count globalFailures =
        globalSum(static_cast<Count>(localFailures));

    return globalFailures == 0 ? 0 : 1;
}

} // namespace

// *********************************** main ***********************************

int main(int argc, char* argv[])
{
    Catch::Session session;

    const int cliResult = session.applyCommandLine(argc, argv);

    if (cliResult != 0)
    {
        return cliResult;
    }

    // Every rank must run the test cases in the same order
    Catch::ConfigData& config = session.configData();
    config.runOrder = Catch::TestRunOrder::Declared;
    config.rngSeed = 1u;

    // Run the session directly without constructing the runtimes
    if (config.listTests || config.listTags || config.listReporters)
    {
        return session.run();
    }

    return runWithRuntimes(session);
}