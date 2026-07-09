/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SubmeshData.h
 * @brief Flat wire format of one rank's submesh
 *
 * @details Everything one rank needs to rebuild its piece of the mesh,
 * as plain vectors: MeshDecomposer fills one SubmeshData per rank on the
 * partitioning rank, MeshDistributor ships it, and MeshCreator turns it
 * back into a Mesh. Variable-length lists (face nodes, cell faces) use
 * CSR layout: offsets[i]..offsets[i+1] slices the flat array.
 *
 * All indices are LOCAL to the submesh. Face order is
 * [internal | processor, grouped by neighbor rank | boundary, by patch],
 * so every patch is a contiguous face range. Cell order is
 * [owned | ghost, grouped by neighbor rank]. A face with no neighbor
 * stores noNeighbor. Geometry is recomputed from nodes on the receiving
 * rank; only ghost cells ship theirs (centroid/volume from the complete
 * mesh), since their remaining faces live on the owning rank.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <limits>
#include <vector>

// Project headers
#include "Integer.h"
#include "Scalar.h"
#include "StringTypes.h"

// ***************************** struct SubmeshData ****************************

struct SubmeshData
{
    /// faceNeighbor value marking a boundary face
    static constexpr Index noNeighbor = std::numeric_limits<Index>::max();

    /// Cells this rank owns
    Count numOwnedCells = 0;

    /// Cell count of the complete mesh (validation reference)
    Count totalCellCount = 0;

    /// Global index of this rank's first owned cell (rank-major)
    Index cellGlobalOffset = 0;

    /// Node coordinates, flat xyz triplets
    ScalarList nodeCoords;

    /// Face -> node connectivity, CSR offsets (numFaces + 1)
    IndexList faceNodeOffsets;

    /// Face -> node connectivity, flat local node indices
    IndexList faceNodes;

    /// Owner cell per face (owned or ghost local index)
    IndexList faceOwner;

    /// Neighbor cell per face, noNeighbor on boundary faces
    IndexList faceNeighbor;

    /// Cell -> face connectivity, CSR offsets (numOwnedCells + 1)
    IndexList cellFaceOffsets;

    /// Cell -> face connectivity, flat local face indices
    IndexList cellFaces;

    /// Face signs aligned with cellFaces (+1 owner side, -1 neighbor)
    std::vector<int8_t> cellFaceSigns;

    /// Physical patch fragments present on this rank
    std::vector<Name> patchNames;

    /// Zone id of each patch fragment
    IndexList patchZoneIds;

    /// PatchType of each patch fragment, stored as its underlying value
    IndexList patchTypes;

    /// First local face of each patch fragment
    IndexList patchFirstFace;

    /// Last local face of each patch fragment
    IndexList patchLastFace;

    /// Neighbor rank of each processor patch, ascending
    IndexList procNeighborRanks;

    /// First local cut face of each processor patch
    IndexList procFirstFace;

    /// Last local cut face of each processor patch
    IndexList procLastFace;

    /// Send cells per processor patch, CSR offsets
    IndexList procSendOffsets;

    /// Send cells, flat owned local indices in the neighbor's ghost order
    IndexList procSendCells;

    /// Ghost cells per processor patch, CSR offsets into the ghost tail
    IndexList procGhostOffsets;

    /// Ghost centroids, flat xyz triplets (from the complete mesh)
    ScalarList ghostCentroids;

    /// Ghost volumes (from the complete mesh)
    ScalarList ghostVolumes;

    /// Global cell index of each ghost
    IndexList ghostGlobalIds;
};
