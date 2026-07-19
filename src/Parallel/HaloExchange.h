/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HaloExchange.h
 * @brief Update ghost cells from their owning ranks
 *
 * @details One exchangeHalos call sends one packed non-blocking message
 * per neighbor rank and fills the ghost tail with the received values.
 * Call it after writing owned cells that neighbors read across a cut.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <cstring>
#include <initializer_list>
#include <vector>

// External library headers
#include <mpi.h>

// Project headers
#include "CellData.h"
#include "Mesh.h"

// ***************************** Halo Exchange ********************************

template<typename T>
void exchangeHalos
(
    const Mesh& mesh,
    std::initializer_list<CellData<T>*> fields
)
{
    constexpr int haloTag = 50;

    const ProcessorPatchList& patches = mesh.processorPatches();

    if (patches.empty())
    {
        return;
    }

    const Count numFields = fields.size();
    const Count numPatches = patches.size();

    std::vector<std::vector<T>> sendBuffers(numPatches);
    std::vector<std::vector<T>> recvBuffers(numPatches);
    std::vector<MPI_Request> requests;
    requests.reserve(2 * numPatches);

    // Receives first, message layout is [field0 cells | field1 cells | ...]
    for (Index p = 0; p < numPatches; ++p)
    {
        recvBuffers[p].resize(numFields * patches[p].ghostCellCount());

        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Irecv
        (
            recvBuffers[p].data(),
            static_cast<int>(recvBuffers[p].size() * sizeof(T)),
            MPI_BYTE,
            static_cast<int>(patches[p].neighborRank()),
            haloTag,
            MPI_COMM_WORLD,
            &request
        );
        requests.push_back(request);
    }

    for (Index p = 0; p < numPatches; ++p)
    {
        const IndexListRef sendCells = patches[p].sendCellIndices();

        sendBuffers[p].reserve(numFields * sendCells.size());

        for (const CellData<T>* field : fields)
        {
            for (const Index cellIdx : sendCells)
            {
                sendBuffers[p].push_back((*field)[cellIdx]);
            }
        }

        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Isend
        (
            sendBuffers[p].data(),
            static_cast<int>(sendBuffers[p].size() * sizeof(T)),
            MPI_BYTE,
            static_cast<int>(patches[p].neighborRank()),
            haloTag,
            MPI_COMM_WORLD,
            &request
        );
        requests.push_back(request);
    }

    MPI_Waitall
    (
        static_cast<int>(requests.size()),
        requests.data(),
        MPI_STATUSES_IGNORE
    );

    // Ghost cells of one patch are contiguous per field: one copy each
    for (Index p = 0; p < numPatches; ++p)
    {
        const Count ghostCount = patches[p].ghostCellCount();
        const Index ghostFirst = patches[p].ghostFirstCell();

        Index f = 0;
        for (CellData<T>* field : fields)
        {
            std::memcpy
            (
                field->data() + ghostFirst,
                recvBuffers[p].data() + f * ghostCount,
                ghostCount * sizeof(T)
            );
            ++f;
        }
    }
}
