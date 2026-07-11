/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HaloExchange.cpp
 * @brief Non-blocking per-neighbor ghost updates
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "HaloExchange.h"

// Standard library headers
#include <cstring>
#include <vector>

// External library headers
#include <mpi.h>

// ****************************** Internal Helpers ****************************

namespace
{

constexpr int haloTag = 50;


template<typename T>
void exchangeFields
(
    const Mesh& mesh,
    const std::vector<CellData<T>*>& fields
)
{
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

    // Post every receive first, then every packed send; message layout
    // is [field0 cells | field1 cells | ...] on both sides
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

        for (Index f = 0; f < numFields; ++f)
        {
            const CellData<T>& field = *fields[f];

            for (const Index cellIdx : sendCells)
            {
                sendBuffers[p].push_back(field[cellIdx]);
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

        for (Index f = 0; f < numFields; ++f)
        {
            std::memcpy
            (
                fields[f]->data() + ghostFirst,
                recvBuffers[p].data() + f * ghostCount,
                ghostCount * sizeof(T)
            );
        }
    }
}

} // namespace

// ***************************** Halo Exchange ********************************

void exchangeHalos
(
    const Mesh& mesh,
    std::initializer_list<ScalarField*> fields
)
{
    exchangeFields(mesh, std::vector<ScalarField*>(fields));
}


void exchangeHalos
(
    const Mesh& mesh,
    VectorField& field
)
{
    exchangeFields(mesh, std::vector<VectorField*>{&field});
}


void exchangeHalos
(
    const Mesh& mesh,
    std::initializer_list<VectorField*> fields
)
{
    exchangeFields(mesh, std::vector<VectorField*>(fields));
}


void exchangeHalos
(
    const Mesh& mesh,
    TensorField& field
)
{
    exchangeFields(mesh, std::vector<TensorField*>{&field});
}
