/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshCreator.cpp
 * @brief Mesh read and geometry-preparation phase, serial or distributed
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "MeshCreator.h"

// Standard library headers
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Project headers
#include "Comm.h"
#include "DecompositionChecker.h"
#include "Logger.h"
#include "MeshChecker.h"
#include "MeshDecomposer.h"
#include "MeshDistributor.h"
#include "MeshReader.h"
#include "SubmeshData.h"

// ****************************** Internal Helpers ****************************

namespace
{

/// Compute face and owned-cell geometry; ghost stubs already carry theirs
void prepareGeometry(Mesh& mesh, bool debug)
{
    std::vector<FaceIntegrals> faceIntegrals(mesh.numFaces());
    auto& faces = mesh.faces();
    const auto& nodes = mesh.nodes();

    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        faceIntegrals[faceIdx] =
            faces[faceIdx].geometricProperties(nodes);
    }
    if (debug)
    {
        std::cout
            << "Geometric properties calculated for faces." << '\n';
    }

    auto& cells = mesh.cells();
    const Count numOwnedCells = mesh.numOwnedCells();

    for (Index cellIdx = 0; cellIdx < numOwnedCells; ++cellIdx)
    {
        cells[cellIdx].geometricProperties(faces, faceIntegrals);
    }
    if (debug)
    {
        std::cout
            << "Geometric properties calculated for cells." << '\n';
    }

    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        faces[faceIdx].distances(cells);
    }
    if (debug)
    {
        std::cout
            << "Distance properties calculated for faces." << '\n';
    }
}


/// Read the complete mesh from file, prepare it, optionally check it
Mesh readCompleteMesh(const CaseConfiguration& config)
{
    MeshReader meshReader(config.meshFile);

    Mesh mesh
    (
        meshReader.moveNodes(),
        meshReader.moveFaces(),
        meshReader.moveCells(),
        meshReader.moveBoundaryPatches()
    );

    std::cout
        << "Mesh Loaded: " << mesh.numNodes() << " nodes, "
        << mesh.numFaces() << " faces, " << mesh.numCells()
        << " cells." << '\n';

    prepareGeometry(mesh, config.debug);

    if (config.checkQuality)
    {
        MeshChecker::check(mesh);
    }

    return mesh;
}


/// Rebuild one rank's Mesh from its flat submesh block
Mesh buildSubmesh(SubmeshData block, bool debug)
{
    const Count numFaces = block.faceOwner.size();
    const Count numGhosts = block.ghostVolumes.size();

    // Nodes
    NodeList nodes;
    nodes.reserve(block.nodeCoords.size() / 3);

    for (Index i = 0; i < block.nodeCoords.size(); i += 3)
    {
        nodes.emplace_back
        (
            block.nodeCoords[i],
            block.nodeCoords[i + 1],
            block.nodeCoords[i + 2]
        );
    }

    // Faces, with original node order and owner/neighbor roles
    FaceList faces;
    faces.reserve(numFaces);

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        IndexList faceNodes
        (
            block.faceNodes.data() + block.faceNodeOffsets[faceIdx],
            block.faceNodes.data() + block.faceNodeOffsets[faceIdx + 1]
        );

        if (block.faceNeighbor[faceIdx] == SubmeshData::noNeighbor)
        {
            faces.emplace_back
            (
                faceIdx,
                std::move(faceNodes),
                block.faceOwner[faceIdx]
            );
        }
        else
        {
            faces.emplace_back
            (
                faceIdx,
                std::move(faceNodes),
                block.faceOwner[faceIdx],
                block.faceNeighbor[faceIdx]
            );
        }
    }

    // Neighbor lists rebuilt from local faces, ghosts included
    std::vector<IndexList> cellNeighbors(block.numOwnedCells);

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        const Index neighbor = block.faceNeighbor[faceIdx];

        if (neighbor == SubmeshData::noNeighbor)
        {
            continue;
        }

        const Index owner = block.faceOwner[faceIdx];

        if (owner < block.numOwnedCells)
        {
            cellNeighbors[owner].push_back(neighbor);
        }
        if (neighbor < block.numOwnedCells)
        {
            cellNeighbors[neighbor].push_back(owner);
        }
    }

    CellList cells;
    cells.reserve(block.numOwnedCells + numGhosts);

    for (Index cellIdx = 0; cellIdx < block.numOwnedCells; ++cellIdx)
    {
        IndexList cellFaces
        (
            block.cellFaces.data() + block.cellFaceOffsets[cellIdx],
            block.cellFaces.data() + block.cellFaceOffsets[cellIdx + 1]
        );
        Cell::FaceSignList signs
        (
            block.cellFaceSigns.data() + block.cellFaceOffsets[cellIdx],
            block.cellFaceSigns.data() + block.cellFaceOffsets[cellIdx + 1]
        );

        cells.emplace_back
        (
            cellIdx,
            std::move(cellFaces),
            std::move(cellNeighbors[cellIdx]),
            std::move(signs)
        );
    }

    // Ghost stubs: geometry from the complete mesh, no local topology
    for (Index i = 0; i < numGhosts; ++i)
    {
        Cell ghost;
        ghost.setIdx(block.numOwnedCells + i);
        ghost.setGeometry
        (
            Vector
            (
                block.ghostCentroids[3 * i],
                block.ghostCentroids[3 * i + 1],
                block.ghostCentroids[3 * i + 2]
            ),
            block.ghostVolumes[i]
        );
        cells.push_back(std::move(ghost));
    }

    // One patch per cut, so the cuts look like ordinary named patches
    PatchList patches;

    for (Index p = 0; p < block.patchNames.size(); ++p)
    {
        BoundaryPatch patch
        (
            block.patchZoneIds[p],
            block.patchFirstFace[p],
            block.patchLastFace[p]
        );
        patch.setName(block.patchNames[p]);
        patch.setType(static_cast<PatchType>(block.patchTypes[p]));
        patches.push_back(std::move(patch));
    }

    ProcessorPatchList processorPatches;

    for (Index p = 0; p < block.procNeighborRanks.size(); ++p)
    {
        BoundaryPatch patch
        (
            0,
            block.procFirstFace[p],
            block.procLastFace[p]
        );
        patch.setName
        (
            "processor" + std::to_string(Comm::myProcessorNum())
          + "to" + std::to_string(block.procNeighborRanks[p])
        );
        patch.setType(PatchType::processor);
        patches.push_back(std::move(patch));

        processorPatches.emplace_back
        (
            block.procNeighborRanks[p],
            block.procFirstFace[p],
            block.procLastFace[p],
            IndexList
            (
                block.procSendCells.data() + block.procSendOffsets[p],
                block.procSendCells.data() + block.procSendOffsets[p + 1]
            ),
            block.numOwnedCells + block.procGhostOffsets[p],
            block.procGhostOffsets[p + 1] - block.procGhostOffsets[p]
        );
    }

    Mesh mesh
    (
        std::move(nodes),
        std::move(faces),
        std::move(cells),
        std::move(patches),
        block.numOwnedCells,
        std::move(block.ghostGlobalIds),
        std::move(processorPatches)
    );

    prepareGeometry(mesh, debug);

    return mesh;
}

} // namespace

// *************************** namespace MeshCreator **************************

namespace MeshCreator
{

Mesh create(const CaseConfiguration& config)
{
    std::cout << '\n';
    Logger::sectionHeader("Reading and Preparing Mesh");

    if (!Comm::parallelRun())
    {
        return readCompleteMesh(config);
    }

    // The master rank reads and partitions, then ships the submeshes
    std::vector<SubmeshData> blocks;

    if (Comm::master())
    {
        {
            const Mesh completeMesh = readCompleteMesh(config);
            const MeshDecomposer decomposer
            (
                completeMesh,
                Comm::numProcessors()
            );

            std::cout
                << "Decomposing into " << Comm::numProcessors()
                << " submeshes (METIS)." << '\n';

            blocks = decomposer.decompose();
        }

        // The complete mesh is gone; the local submesh takes its place
        Mesh::resetCounts();
    }

    SubmeshData block = MeshDistributor::distribute(std::move(blocks));

    const Count totalCellCount = block.totalCellCount;

    Mesh mesh = buildSubmesh(std::move(block), config.debug);

    DecompositionChecker::check(mesh, totalCellCount);

    return mesh;
}

} // namespace MeshCreator
