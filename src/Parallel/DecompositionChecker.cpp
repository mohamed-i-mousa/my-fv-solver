/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file DecompositionChecker.cpp
 * @brief Exchange-and-compare validation across every inter-rank cut
 *****************************************************************************/

// ********************************** Headers *********************************

#include "DecompositionChecker.h"

// Standard library headers
#include <cmath>
#include <string>
#include <vector>

// External library headers
#include <mpi.h>

// Project headers
#include "BoundaryConditions.h"
#include "Comm.h"
#include "ErrorHandler.h"
#include "Field.h"
#include "LeastSquares.h"
#include "Logger.h"
#include "MpiScalarType.h"
#include "Reduce.h"

// ****************************** Internal Helpers ****************************

namespace
{

constexpr int exchangeTag = 42;


/// Swap scalar payloads with the neighbor across one cut; the two sides
/// may send different amounts (my send count equals the NEIGHBOR's ghost
/// count, not my own — one owned cell can face several remote cells)
std::vector<Scalar> exchangeWithNeighbor
(
    const ProcessorPatch& patch,
    const std::vector<Scalar>& sendBuffer,
    Count recvCount
)
{
    std::vector<Scalar> recvBuffer(recvCount);

    // Sendrecv is cycle-safe: send and receive progress concurrently, so
    // per-rank patch ordering cannot deadlock the neighbor loop
    MPI_Sendrecv
    (
        sendBuffer.data(),
        static_cast<int>(sendBuffer.size()),
        mpiScalarType(),
        static_cast<int>(patch.neighborRank()),
        exchangeTag,
        recvBuffer.data(),
        static_cast<int>(recvBuffer.size()),
        mpiScalarType(),
        static_cast<int>(patch.neighborRank()),
        exchangeTag,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    return recvBuffer;
}


/// Check the owned-cell counts sum to the complete mesh
void checkCellAccounting
(
    const Mesh& mesh,
    Count totalCellCount
)
{
    const Count ownedTotal = globalSum(mesh.numOwnedCells());

    if (ownedTotal != totalCellCount)
    {
        FatalError
        (
            "Decomposition check: owned cells sum to "
          + std::to_string(ownedTotal) + ", complete mesh has "
          + std::to_string(totalCellCount)
        );
    }
}


/// Ghost centroids/volumes must BIT-match the owning rank's cells
void checkGhostGeometry(const Mesh& mesh)
{
    const CellListRef cells = mesh.cells();

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        // My owned send cells feed the neighbor's ghosts; the neighbor's
        // sends arrive in exactly my ghost order (the ordering contract)
        std::vector<Scalar> sendBuffer;
        sendBuffer.reserve(4 * patch.sendCellIndices().size());

        for (const Index cellIdx : patch.sendCellIndices())
        {
            const Vector& centroid = cells[cellIdx].centroid();
            sendBuffer.push_back(centroid.x());
            sendBuffer.push_back(centroid.y());
            sendBuffer.push_back(centroid.z());
            sendBuffer.push_back(cells[cellIdx].volume());
        }

        const std::vector<Scalar> received =
            exchangeWithNeighbor
            (
                patch,
                sendBuffer,
                4 * patch.ghostCellCount()
            );

        for (Index i = 0; i < patch.ghostCellCount(); ++i)
        {
            const Cell& ghost = cells[patch.ghostFirstCell() + i];

            const bool match =
                received[4 * i]     == ghost.centroid().x()
             && received[4 * i + 1] == ghost.centroid().y()
             && received[4 * i + 2] == ghost.centroid().z()
             && received[4 * i + 3] == ghost.volume();

            if (!match)
            {
                FatalError
                (
                    "Decomposition check: ghost cell "
                  + std::to_string(patch.ghostFirstCell() + i)
                  + " does not bit-match its owner on rank "
                  + std::to_string(patch.neighborRank())
                );
            }
        }
    }
}


/// Both copies of every cut face must BIT-match in geometry
void checkCutFaceGeometry(const Mesh& mesh)
{
    const FaceListRef faces = mesh.faces();

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        std::vector<Scalar> sendBuffer;
        sendBuffer.reserve(7 * patch.numFaces());

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            const Face& face = faces[faceIdx];
            sendBuffer.push_back(face.centroid().x());
            sendBuffer.push_back(face.centroid().y());
            sendBuffer.push_back(face.centroid().z());
            sendBuffer.push_back(face.normal().x());
            sendBuffer.push_back(face.normal().y());
            sendBuffer.push_back(face.normal().z());
            sendBuffer.push_back(face.projectedArea());
        }

        const std::vector<Scalar> received =
            exchangeWithNeighbor(patch, sendBuffer, sendBuffer.size());

        for (Index i = 0; i < sendBuffer.size(); ++i)
        {
            if (received[i] != sendBuffer[i])
            {
                FatalError
                (
                    "Decomposition check: cut face geometry differs from "
                    "rank " + std::to_string(patch.neighborRank())
                  + " (both sides recompute from identical node data, so "
                    "any difference is an extraction bug)"
                );
            }
        }
    }
}


/// Least-squares gradient of a linear field, exact across the cuts
void checkGradientAcrossCuts(const Mesh& mesh)
{
    // phi = a x + b y + c z + d sampled at every centroid, ghosts
    // included — the ghost values stand in for the halo exchange
    const Scalar a = S(1.5);
    const Scalar b = S(2.5);
    const Scalar c = S(-3.5);
    const Scalar d = S(4.2);

    ScalarField phi;

    const CellListRef cells = mesh.cells();
    const Count numCells = mesh.numCells();

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Vector& x = cells[cellIdx].centroid();
        phi[cellIdx] = a * x.x() + b * x.y() + c * x.z() + d;
    }

    // Geometric-only construction: no boundary values are ever read,
    // because cells touching a physical boundary are skipped below
    const BoundaryConditions bc;
    const LeastSquares leastSquares(mesh, bc);

    // Owned cells adjacent to a cut, away from physical boundaries
    std::vector<bool> tested(mesh.numOwnedCells(), false);
    const FaceListRef faces = mesh.faces();

    Scalar maxError = S(0.0);
    Count testedCells = 0;

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            const Face& face = faces[faceIdx];
            const Index neighborIdx = face.neighborCell().value();

            const Index ownedSide =
                face.ownerCell() < mesh.numOwnedCells()
              ? face.ownerCell()
              : neighborIdx;

            if (tested[ownedSide])
            {
                continue;
            }

            tested[ownedSide] = true;

            bool touchesBoundary = false;

            for (const Index cellFaceIdx : cells[ownedSide].faceIndices())
            {
                if (faces[cellFaceIdx].isBoundary())
                {
                    touchesBoundary = true;
                    break;
                }
            }

            if (touchesBoundary)
            {
                continue;
            }

            const Vector gradient =
                leastSquares.cellGradient(Field::p, phi, ownedSide);

            maxError =
                std::max
                (
                    maxError,
                    std::max
                    (
                        std::abs(gradient.x() - a),
                        std::max
                        (
                            std::abs(gradient.y() - b),
                            std::abs(gradient.z() - c)
                        )
                    )
                );

            ++testedCells;
        }
    }

    const Scalar globalError = globalMax(maxError);
    const Count globalTested = globalSum(testedCells);

    // A wrong or stale ghost stencil produces O(1) errors; a healthy
    // least-squares reconstruction of a linear field is exact to
    // round-off amplified by the stencil conditioning
    if (globalError > S(1e-8))
    {
        FatalError
        (
            "Decomposition check: least-squares gradient across the cuts "
            "reached error " + std::to_string(globalError)
        );
    }

    if (Comm::master())
    {
        Logger::keyValue("Gradient test cells", globalTested);
        Logger::keyValue("Gradient max error", globalError);
    }
}

} // namespace

// ********************** namespace DecompositionChecker **********************

void DecompositionChecker::check
(
    const Mesh& mesh,
    Count totalCellCount
)
{
    checkCellAccounting(mesh, totalCellCount);
    checkGhostGeometry(mesh);
    checkCutFaceGeometry(mesh);

    Count numCutFaces = 0;

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        numCutFaces += patch.numFaces();
    }

    const Count globalGhosts = globalSum(mesh.numGhostCells());
    const Count globalCuts = globalSum(numCutFaces);

    if (Comm::master())
    {
        Logger::subsection("Decomposition checks");
        Logger::keyValue("Total owned cells", totalCellCount);
        Logger::keyValue("Ghost cells (all ranks)", globalGhosts);
        Logger::keyValue("Cut faces (both sides)", globalCuts);
        Logger::keyValue("Ghost geometry", "bit-identical");
        Logger::keyValue("Cut-face geometry", "bit-identical");
    }

    checkGradientAcrossCuts(mesh);
}
