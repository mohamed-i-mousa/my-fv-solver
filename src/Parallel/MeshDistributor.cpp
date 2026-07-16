/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshDistributor.cpp
 * @brief Byte-stream serialization and MPI scatter of submesh blocks
 *****************************************************************************/

// ********************************** Headers *********************************

#include "MeshDistributor.h"

// Standard library headers
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>

// External library headers
#include <mpi.h>

// Project headers
#include "Comm.h"
#include "ErrorHandler.h"

// ****************************** Internal Helpers ****************************

namespace
{

// MPI message tags for the two-message (length, payload) handshake
constexpr int sizeTag = 40;
constexpr int payloadTag = 41;

/// Append-only byte buffer with typed writes
class ByteWriter
{
public:

    template<typename T>
    void write(T value)
    {
        static_assert(std::is_trivially_copyable_v<T>);

        const Count at = buffer_.size();
        buffer_.resize(at + sizeof(T));
        std::memcpy(buffer_.data() + at, &value, sizeof(T));
    }

    template<typename T>
    void writeList(const std::vector<T>& values)
    {
        static_assert(std::is_trivially_copyable_v<T>);

        write<Count>(values.size());

        const Count at = buffer_.size();
        buffer_.resize(at + values.size() * sizeof(T));
        std::memcpy
        (
            buffer_.data() + at,
            values.data(),
            values.size() * sizeof(T)
        );
    }

    void writeNames(const std::vector<Name>& names)
    {
        write<Count>(names.size());

        for (const Name& name : names)
        {
            write<Count>(name.size());

            const Count at = buffer_.size();
            buffer_.resize(at + name.size());
            std::memcpy(buffer_.data() + at, name.data(), name.size());
        }
    }

    [[nodiscard]] const std::vector<char>& buffer() const noexcept
    {
        return buffer_;
    }

private:

    std::vector<char> buffer_;
};


/// Sequential typed reads over a received byte buffer
class ByteReader
{
public:

    explicit ByteReader(std::vector<char> buffer)
    :
        buffer_(std::move(buffer))
    {}

    template<typename T>
    [[nodiscard]] T read()
    {
        static_assert(std::is_trivially_copyable_v<T>);

        checkRemaining(sizeof(T));

        T value;
        std::memcpy(&value, buffer_.data() + cursor_, sizeof(T));
        cursor_ += sizeof(T);

        return value;
    }

    template<typename T>
    [[nodiscard]] std::vector<T> readList()
    {
        static_assert(std::is_trivially_copyable_v<T>);

        const Count size = read<Count>();
        checkRemaining(size * sizeof(T));

        std::vector<T> values(size);
        std::memcpy(values.data(), buffer_.data() + cursor_,
                    size * sizeof(T));
        cursor_ += size * sizeof(T);

        return values;
    }

    [[nodiscard]] std::vector<Name> readNames()
    {
        const Count size = read<Count>();

        std::vector<Name> names;
        names.reserve(size);

        for (Index i = 0; i < size; ++i)
        {
            const Count length = read<Count>();
            checkRemaining(length);

            names.emplace_back(buffer_.data() + cursor_, length);
            cursor_ += length;
        }

        return names;
    }

private:

    void checkRemaining(Count bytes) const
    {
        if (cursor_ + bytes > buffer_.size())
        {
            FatalError
            (
                "MeshDistributor: submesh payload shorter than its layout"
            );
        }
    }

    std::vector<char> buffer_;
    Count cursor_ = 0;
};


/// Flatten one block; field order is the deserialize() contract
std::vector<char> serialize(const SubmeshData& block)
{
    ByteWriter writer;

    writer.write<Count>(block.numOwnedCells);
    writer.write<Count>(block.totalCellCount);
    writer.write<Index>(block.cellGlobalOffset);
    writer.writeList(block.nodeCoords);
    writer.writeList(block.faceNodeOffsets);
    writer.writeList(block.faceNodes);
    writer.writeList(block.faceOwner);
    writer.writeList(block.faceNeighbor);
    writer.writeList(block.cellFaceOffsets);
    writer.writeList(block.cellFaces);
    writer.writeList(block.cellFaceSigns);
    writer.writeNames(block.patchNames);
    writer.writeList(block.patchZoneIds);
    writer.writeList(block.patchTypes);
    writer.writeList(block.patchFirstFace);
    writer.writeList(block.patchLastFace);
    writer.writeList(block.procNeighborRanks);
    writer.writeList(block.procFirstFace);
    writer.writeList(block.procLastFace);
    writer.writeList(block.procSendOffsets);
    writer.writeList(block.procSendCells);
    writer.writeList(block.procGhostOffsets);
    writer.writeList(block.ghostCentroids);
    writer.writeList(block.ghostVolumes);
    writer.writeList(block.ghostGlobalIds);

    return writer.buffer();
}


/// Rebuild one block; mirrors serialize() field for field
SubmeshData deserialize(std::vector<char> bytes)
{
    ByteReader reader(std::move(bytes));

    SubmeshData block;

    block.numOwnedCells = reader.read<Count>();
    block.totalCellCount = reader.read<Count>();
    block.cellGlobalOffset = reader.read<Index>();
    block.nodeCoords = reader.readList<Scalar>();
    block.faceNodeOffsets = reader.readList<Index>();
    block.faceNodes = reader.readList<Index>();
    block.faceOwner = reader.readList<Index>();
    block.faceNeighbor = reader.readList<Index>();
    block.cellFaceOffsets = reader.readList<Index>();
    block.cellFaces = reader.readList<Index>();
    block.cellFaceSigns = reader.readList<int8_t>();
    block.patchNames = reader.readNames();
    block.patchZoneIds = reader.readList<Index>();
    block.patchTypes = reader.readList<Index>();
    block.patchFirstFace = reader.readList<Index>();
    block.patchLastFace = reader.readList<Index>();
    block.procNeighborRanks = reader.readList<Index>();
    block.procFirstFace = reader.readList<Index>();
    block.procLastFace = reader.readList<Index>();
    block.procSendOffsets = reader.readList<Index>();
    block.procSendCells = reader.readList<Index>();
    block.procGhostOffsets = reader.readList<Index>();
    block.ghostCentroids = reader.readList<Scalar>();
    block.ghostVolumes = reader.readList<Scalar>();
    block.ghostGlobalIds = reader.readList<Index>();

    return block;
}

} // namespace

// ************************* namespace MeshDistributor ************************

SubmeshData MeshDistributor::distribute(std::vector<SubmeshData> blocks)
{
    if (Comm::master())
    {
        if (blocks.size() != Comm::numProcessors())
        {
            FatalError
            (
                "MeshDistributor: expected one submesh block per rank"
            );
        }

        for (Index rank = 1; rank < Comm::numProcessors(); ++rank)
        {
            const std::vector<char> bytes = serialize(blocks[rank]);

            if (bytes.size()
              > static_cast<Count>(std::numeric_limits<int>::max()))
            {
                FatalError
                (
                    "MeshDistributor: submesh block exceeds the 2 GiB "
                    "single-message limit"
                );
            }

            const unsigned long long size = bytes.size();
            MPI_Send
            (
                &size,
                1,
                MPI_UNSIGNED_LONG_LONG,
                static_cast<int>(rank),
                sizeTag,
                MPI_COMM_WORLD
            );
            MPI_Send
            (
                bytes.data(),
                static_cast<int>(bytes.size()),
                MPI_BYTE,
                static_cast<int>(rank),
                payloadTag,
                MPI_COMM_WORLD
            );
        }

        return std::move(blocks[0]);
    }

    unsigned long long size = 0;
    MPI_Recv
    (
        &size,
        1,
        MPI_UNSIGNED_LONG_LONG,
        0,
        sizeTag,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    std::vector<char> bytes(size);
    MPI_Recv
    (
        bytes.data(),
        static_cast<int>(size),
        MPI_BYTE,
        0,
        payloadTag,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    return deserialize(std::move(bytes));
}
