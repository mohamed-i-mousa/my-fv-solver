/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshDecomposer.cpp
 * @brief METIS partitioning and per-rank submesh extraction
 *****************************************************************************/

// ********************************** Headers *********************************

#include "MeshDecomposer.h"

// Standard library headers
#include <limits>
#include <string>

// External library headers
#include <metis.h>

// Project headers
#include "ErrorHandler.h"

// **************************** Build-Time Invariants **************************

// The dual graph is handed to METIS as idx_t arrays; the 32-bit build is
// assumed throughout (matches PetscInt and the < 2^31 cells decision)
static_assert
(
    IDXTYPEWIDTH == 32,
    "MeshDecomposer assumes a 32-bit METIS idx_t build"
);

// ****************************** Internal Helpers ****************************

namespace
{

/// Sentinel for "not present on this rank" index maps
constexpr Index invalidIdx = std::numeric_limits<Index>::max();

} // namespace

// ************************* Special Member Functions *************************

MeshDecomposer::MeshDecomposer
(
    const Mesh& mesh,
    Count numParts
)
:
    mesh_(mesh),
    numParts_(numParts)
{
    if (numParts_ == 0)
    {
        FatalError("MeshDecomposer: numParts must be at least 1");
    }

    if (mesh_.numCells() >= std::numeric_limits<idx_t>::max())
    {
        FatalError
        (
            "MeshDecomposer: cell count exceeds the 32-bit index range"
        );
    }
}

// ****************************** Public Methods ******************************

IndexList MeshDecomposer::partition() const
{
    const Count numCells = mesh_.numCells();

    IndexList cellToRank(numCells, 0);

    if (numParts_ == 1)
    {
        return cellToRank;
    }

    // Dual graph in CSR form: cells are vertices, internal faces edges
    const FaceListRef faces = mesh_.faces();

    std::vector<idx_t> xadj(numCells + 1, 0);

    Count numInternalFaces = 0;

    for (const Face& face : faces)
    {
        if (!face.isBoundary())
        {
            ++xadj[face.ownerCell() + 1];
            ++xadj[face.neighborCell().value() + 1];
            ++numInternalFaces;
        }
    }

    // Two CSR entries per internal face, overflowing idx_t before cells do
    if 
    (
        2 * numInternalFaces
     >= static_cast<Count>(std::numeric_limits<idx_t>::max())
    )
    {
        FatalError
        (
            "MeshDecomposer: adjacency entry count exceeds the 32-bit "
            "index range"
        );
    }

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        xadj[cellIdx + 1] += xadj[cellIdx];
    }

    std::vector<idx_t> adjncy(static_cast<Count>(xadj[numCells]));
    std::vector<idx_t> cursor(xadj.begin(), xadj.end() - 1);

    for (const Face& face : faces)
    {
        if (!face.isBoundary())
        {
            const Index owner = face.ownerCell();
            const Index neighbor = face.neighborCell().value();

            adjncy[static_cast<Count>(cursor[owner]++)] =
                static_cast<idx_t>(neighbor);
            adjncy[static_cast<Count>(cursor[neighbor]++)] =
                static_cast<idx_t>(owner);
        }
    }

    idx_t numVertices = static_cast<idx_t>(numCells);
    idx_t numConstraints = 1;
    idx_t numParts = static_cast<idx_t>(numParts_);
    idx_t edgeCut = 0;
    std::vector<idx_t> part(numCells, 0);

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // Fixed seed: every run (and every rank, if ever run redundantly)
    // must produce the identical partition
    options[METIS_OPTION_SEED] = 0;
    options[METIS_OPTION_NUMBERING] = 0;

    const int status =
        METIS_PartGraphKway
        (
            &numVertices,
            &numConstraints,
            xadj.data(),
            adjncy.data(),
            nullptr,            // vertex weights: uniform
            nullptr,            // vertex sizes: unused
            nullptr,            // edge weights: uniform
            &numParts,
            nullptr,            // target part weights: uniform
            nullptr,            // imbalance tolerance: default
            options,
            &edgeCut,
            part.data()
        );

    if (status != METIS_OK)
    {
        FatalError
        (
            "MeshDecomposer: METIS_PartGraphKway failed with status "
          + std::to_string(status)
        );
    }

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        cellToRank[cellIdx] = static_cast<Index>(part[cellIdx]);
    }

    return cellToRank;
}


std::vector<SubmeshData> MeshDecomposer::decompose() const
{
    const Count numCells = mesh_.numCells();
    const Count numFaces = mesh_.numFaces();

    const IndexList cellToRank = partition();

    // Local index of every cell within its rank (rank-major numbering:
    // original-mesh order preserved within each rank)
    IndexList localCellIdx(numCells, 0);
    IndexList rankCellCounts(numParts_, 0);

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        localCellIdx[cellIdx] = rankCellCounts[cellToRank[cellIdx]]++;
    }

    // Rank-major global offsets: rank r's cells follow ranks 0..r-1
    IndexList rankOffsets(numParts_ + 1, 0);

    for (Index rank = 0; rank < numParts_; ++rank)
    {
        rankOffsets[rank + 1] = rankOffsets[rank] + rankCellCounts[rank];
    }

    // Patch of every boundary face (patches are contiguous face ranges)
    IndexList faceToPatch(numFaces, invalidIdx);

    const PatchListRef patches = mesh_.patches();

    for (Index patchIdx = 0; patchIdx < patches.size(); ++patchIdx)
    {
        const BoundaryPatch& patch = patches[patchIdx];

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            faceToPatch[faceIdx] = patchIdx;
        }
    }

    std::vector<SubmeshData> blocks;
    blocks.reserve(numParts_);

    for (Index rank = 0; rank < numParts_; ++rank)
    {
        blocks.push_back
        (
            extractSubmesh
            (
                rank,
                cellToRank,
                localCellIdx,
                rankOffsets,
                faceToPatch
            )
        );
    }

    return blocks;
}

// ****************************** Private Methods *****************************

SubmeshData MeshDecomposer::extractSubmesh
(
    Index rank,
    const IndexList& cellToRank,
    const IndexList& localCellIdx,
    const IndexList& rankOffsets,
    const IndexList& faceToPatch
) const
{
    const FaceListRef faces = mesh_.faces();
    const CellListRef cells = mesh_.cells();
    const PatchListRef patches = mesh_.patches();

    const Count numOwned = rankOffsets[rank + 1] - rankOffsets[rank];

    // ------------------------------------------------------------------
    // Classify every original face touched by this rank, preserving
    // original-mesh order within each group (the ordering contract both
    // sides of every cut derive the send/ghost sequence from)
    // ------------------------------------------------------------------
    IndexList internalFaces;
    std::vector<IndexList> cutFacesTo(numParts_);
    std::vector<IndexList> boundaryFacesOf(patches.size());

    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        const Face& face = faces[faceIdx];
        const bool ownerHere = cellToRank[face.ownerCell()] == rank;

        if (face.isBoundary())
        {
            if (ownerHere)
            {
                if (faceToPatch[faceIdx] == invalidIdx)
                {
                    FatalError
                    (
                        "MeshDecomposer: boundary face "
                      + std::to_string(faceIdx)
                      + " belongs to no patch"
                    );
                }

                boundaryFacesOf[faceToPatch[faceIdx]].push_back(faceIdx);
            }

            continue;
        }

        const bool neighborHere =
            cellToRank[face.neighborCell().value()] == rank;

        if (ownerHere && neighborHere)
        {
            internalFaces.push_back(faceIdx);
        }
        else if (ownerHere || neighborHere)
        {
            const Index remoteCell =
                ownerHere ? face.neighborCell().value() : face.ownerCell();

            cutFacesTo[cellToRank[remoteCell]].push_back(faceIdx);
        }
    }

    // ------------------------------------------------------------------
    // Ghost layer and send lists: walk each cut's faces in order; first
    // appearance fixes both the ghost order here and the send order on
    // the neighbor (which walks the same faces in the same order)
    // ------------------------------------------------------------------
    SubmeshData block;
    block.numOwnedCells = numOwned;
    block.totalCellCount = mesh_.numCells();
    block.cellGlobalOffset = rankOffsets[rank];

    IndexList ghostOf(mesh_.numCells(), invalidIdx);
    IndexList sentTo(mesh_.numCells(), invalidIdx);

    Count numGhosts = 0;

    block.procGhostOffsets.push_back(0);
    block.procSendOffsets.push_back(0);

    for (Index nbr = 0; nbr < numParts_; ++nbr)
    {
        if (cutFacesTo[nbr].empty())
        {
            continue;
        }

        block.procNeighborRanks.push_back(nbr);

        for (const Index faceIdx : cutFacesTo[nbr])
        {
            const Face& face = faces[faceIdx];
            const bool ownerHere = cellToRank[face.ownerCell()] == rank;

            const Index ownCell =
                ownerHere ? face.ownerCell() : face.neighborCell().value();
            const Index remoteCell =
                ownerHere ? face.neighborCell().value() : face.ownerCell();

            if (ghostOf[remoteCell] == invalidIdx)
            {
                ghostOf[remoteCell] = numOwned + numGhosts;
                ++numGhosts;

                const Vector& centroid = cells[remoteCell].centroid();
                block.ghostCentroids.push_back(centroid.x());
                block.ghostCentroids.push_back(centroid.y());
                block.ghostCentroids.push_back(centroid.z());
                block.ghostVolumes.push_back(cells[remoteCell].volume());
                block.ghostGlobalIds.push_back
                (
                    rankOffsets[nbr] + localCellIdx[remoteCell]
                );
            }

            // sentTo marks the neighbor a cell was last queued for, so a
            // cell feeding two cuts is sent to both, but once each
            if (sentTo[ownCell] != nbr)
            {
                sentTo[ownCell] = nbr;
                block.procSendCells.push_back(localCellIdx[ownCell]);
            }
        }

        block.procGhostOffsets.push_back(numGhosts);
        block.procSendOffsets.push_back(block.procSendCells.size());
    }

    // ------------------------------------------------------------------
    // Local face numbering: [internal | cuts by neighbor | boundary by
    // patch] — every patch ends up a contiguous local face range
    // ------------------------------------------------------------------
    IndexList localFaceIdx(faces.size(), invalidIdx);
    IndexList includedFaces;

    for (const Index faceIdx : internalFaces)
    {
        localFaceIdx[faceIdx] = includedFaces.size();
        includedFaces.push_back(faceIdx);
    }

    for (Index patchPos = 0; patchPos < block.procNeighborRanks.size();
         ++patchPos)
    {
        const Index nbr = block.procNeighborRanks[patchPos];

        block.procFirstFace.push_back(includedFaces.size());

        for (const Index faceIdx : cutFacesTo[nbr])
        {
            localFaceIdx[faceIdx] = includedFaces.size();
            includedFaces.push_back(faceIdx);
        }

        block.procLastFace.push_back(includedFaces.size() - 1);
    }

    for (Index patchIdx = 0; patchIdx < patches.size(); ++patchIdx)
    {
        if (boundaryFacesOf[patchIdx].empty())
        {
            continue;
        }

        const BoundaryPatch& patch = patches[patchIdx];

        block.patchNames.push_back(patch.patchName());
        block.patchZoneIds.push_back(patch.zoneIdx());
        block.patchTypes.push_back(static_cast<Index>(patch.type()));
        block.patchFirstFace.push_back(includedFaces.size());

        for (const Index faceIdx : boundaryFacesOf[patchIdx])
        {
            localFaceIdx[faceIdx] = includedFaces.size();
            includedFaces.push_back(faceIdx);
        }

        block.patchLastFace.push_back(includedFaces.size() - 1);
    }

    // ------------------------------------------------------------------
    // Node subset: first appearance along the local face order
    // ------------------------------------------------------------------
    const NodeListRef nodes = mesh_.nodes();

    IndexList localNodeIdx(nodes.size(), invalidIdx);
    Count numLocalNodes = 0;

    for (const Index faceIdx : includedFaces)
    {
        for (const Index nodeIdx : faces[faceIdx].nodeIndices())
        {
            if (localNodeIdx[nodeIdx] == invalidIdx)
            {
                localNodeIdx[nodeIdx] = numLocalNodes++;

                const Vector& coords = nodes[nodeIdx];
                block.nodeCoords.push_back(coords.x());
                block.nodeCoords.push_back(coords.y());
                block.nodeCoords.push_back(coords.z());
            }
        }
    }

    // ------------------------------------------------------------------
    // Faces: original node order and owner/neighbor roles are preserved
    // verbatim (a remote cell simply maps to its ghost index), so both
    // ranks' copies of a cut face carry identical data
    // ------------------------------------------------------------------
    block.faceNodeOffsets.push_back(0);

    for (const Index faceIdx : includedFaces)
    {
        const Face& face = faces[faceIdx];

        for (const Index nodeIdx : face.nodeIndices())
        {
            block.faceNodes.push_back(localNodeIdx[nodeIdx]);
        }

        block.faceNodeOffsets.push_back(block.faceNodes.size());

        const Index owner = face.ownerCell();
        block.faceOwner.push_back
        (
            cellToRank[owner] == rank ? localCellIdx[owner] : ghostOf[owner]
        );

        if (face.isBoundary())
        {
            block.faceNeighbor.push_back(SubmeshData::noNeighbor);
        }
        else
        {
            const Index neighbor = face.neighborCell().value();
            block.faceNeighbor.push_back
            (
                cellToRank[neighbor] == rank
              ? localCellIdx[neighbor]
              : ghostOf[neighbor]
            );
        }
    }

    // ------------------------------------------------------------------
    // Owned cells in local order: face lists and signs copied verbatim
    // (roles unchanged, so the +1 owner / -1 neighbor signs still hold)
    // ------------------------------------------------------------------
    block.cellFaceOffsets.push_back(0);

    IndexList ownedOriginal(numOwned, 0);

    for (Index cellIdx = 0; cellIdx < mesh_.numCells(); ++cellIdx)
    {
        if (cellToRank[cellIdx] == rank)
        {
            ownedOriginal[localCellIdx[cellIdx]] = cellIdx;
        }
    }

    for (const Index cellIdx : ownedOriginal)
    {
        const Cell& cell = cells[cellIdx];
        const IndexListRef cellFaces = cell.faceIndices();
        const Cell::FaceSignListRef signs = cell.faceSigns();

        for (Index j = 0; j < cellFaces.size(); ++j)
        {
            if (localFaceIdx[cellFaces[j]] == invalidIdx)
            {
                FatalError
                (
                    "MeshDecomposer: face of owned cell "
                  + std::to_string(cellIdx)
                  + " missing from its rank's submesh"
                );
            }

            block.cellFaces.push_back(localFaceIdx[cellFaces[j]]);
            block.cellFaceSigns.push_back(signs[j]);
        }

        block.cellFaceOffsets.push_back(block.cellFaces.size());
    }

    return block;
}
