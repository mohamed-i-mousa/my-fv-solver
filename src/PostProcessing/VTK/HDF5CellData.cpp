/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HDF5CellData.cpp
 * @brief Implementation of VTKHDF cell data writer
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "HDF5CellData.h"

// Standard library headers
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

// Project headers
#include "Comm.h"
#include "ErrorHandler.h"
#include "HDF5Wrapper.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

// ***************************** Internal Helpers *****************************

namespace
{

// VTK cell type id for arbitrary polyhedra (VTK_POLYHEDRON)
constexpr unsigned char polyhedronCellType = 42;


Vector newellNormal
(
    IndexListRef faceNodes,
    NodeListRef allNodes
)
{
    Vector normal;

    for (Index i = 0; i < faceNodes.size(); ++i)
    {
        const Vector& current = allNodes[faceNodes[i]];
        const Vector& next = allNodes[faceNodes[(i + 1) % faceNodes.size()]];

        normal += Vector
        (
            (current.y() - next.y()) * (current.z() + next.z()),
            (current.z() - next.z()) * (current.x() + next.x()),
            (current.x() - next.x()) * (current.y() + next.y())
        );
    }

    return normal;
}

} // namespace

// ************************* Special Member Functions *************************

HDF5CellData::HDF5CellData
(
    FilePath fileName,
    const Mesh& mesh,
    bool debug
)
:
    mesh_{mesh},
    fileName_{std::move(fileName)},
    debug_{debug},
    file_{H5I_INVALID_HID},
    vtkhdfGroup_{H5I_INVALID_HID},
    cellDataGroup_{H5I_INVALID_HID},
    stepsGroup_{H5I_INVALID_HID},
    cellDataOffsetsGroup_{H5I_INVALID_HID},
    geometryWritten_{false},
    finalized_{false},
    stepCount_{0},
    numCells_{0},
    numNodes_{0},
    globalNumCells_{0},
    cellRowOffset_{0}
{
    file_ = HDF5::createFile(fileName_);
    vtkhdfGroup_ = HDF5::createGroup(file_, "VTKHDF");

    const long long version[2] = {2, 4};
    HDF5::writeIntArrayAttribute(vtkhdfGroup_, "Version", version, 2);
    HDF5::writeStringAttribute(vtkhdfGroup_, "Type", "UnstructuredGrid");

    cellDataGroup_ = HDF5::createGroup(vtkhdfGroup_, "CellData");
    stepsGroup_ = HDF5::createGroup(vtkhdfGroup_, "Steps");
    cellDataOffsetsGroup_ = HDF5::createGroup(stepsGroup_, "CellDataOffsets");

    HDF5::writeIntScalarAttribute(stepsGroup_, "NSteps", 0);

    // Keep the file closed between writes so it stays readable mid-run
    closeHandles();
}


HDF5CellData::~HDF5CellData() noexcept
{
    finalize();
}

// ****************************** Public Methods ******************************

void HDF5CellData::writeGeometry()
{
    if (geometryWritten_)
    {
        FatalError
        (
            "VTKHDF geometry has already been written: " + fileName_
        );
    }

    openHandles();

    const NodeListRef allNodes = mesh_.nodes();
    const CellListRef allCells = mesh_.cells();
    const FaceListRef allFaces = mesh_.faces();

    // Owned cells only: the neighbor ranks write their own
    const Count numOwnedCells = mesh_.numOwnedCells();

    // Point coordinates, written once as double precision
    std::vector<double> points;
    points.reserve(allNodes.size() * 3);

    for (const Vector& node : allNodes)
    {
        points.push_back(static_cast<double>(node.x()));
        points.push_back(static_cast<double>(node.y()));
        points.push_back(static_cast<double>(node.z()));
    }

    // Per-cell unique point lists plus per-cell-occurrence face streams
    std::vector<long long> connectivity;
    std::vector<long long> offsets;
    const std::vector<unsigned char> types
    (
        numOwnedCells,
        polyhedronCellType
    );
    std::vector<long long> faceConnectivity;
    std::vector<long long> faceOffsets;
    std::vector<long long> polyhedronToFaces;
    std::vector<long long> polyhedronOffsets;

    offsets.reserve(numOwnedCells + 1);
    polyhedronOffsets.reserve(numOwnedCells + 1);
    offsets.push_back(0);
    faceOffsets.push_back(0);
    polyhedronOffsets.push_back(0);

    std::vector<long long> uniquePointIds;
    std::vector<long long> orientedFaceNodes;

    for (Index cellIdx = 0; cellIdx < numOwnedCells; ++cellIdx)
    {
        const Cell& cell = allCells[cellIdx];
        const auto faceIndices = cell.faceIndices();

        uniquePointIds.clear();

        for (Index faceIdx : faceIndices)
        {
            const Face& face = allFaces[faceIdx];
            const auto nodeIndices = face.nodeIndices();

            if (nodeIndices.size() < 3)
            {
                FatalError
                (
                    "Cannot write VTKHDF polyhedron cell "
                  + std::to_string(cellIdx) + ": face "
                  + std::to_string(faceIdx)
                  + " has fewer than three nodes."
                );
            }

            for (Index nodeIdx : nodeIndices)
            {
                const long long pointId = static_cast<long long>(nodeIdx);

                if
                (
                    std::find
                    (
                        uniquePointIds.begin(),
                        uniquePointIds.end(),
                        pointId
                    ) == uniquePointIds.end()
                )
                {
                    uniquePointIds.push_back(pointId);
                }
            }

            orientedFaceNodes.clear();

            for (Index nodeIdx : nodeIndices)
            {
                orientedFaceNodes.push_back(static_cast<long long>(nodeIdx));
            }

            const Vector faceNormal = newellNormal(nodeIndices, allNodes);
            const Vector outwardDirection =
                face.centroid() - cell.centroid();

            if (dot(faceNormal, outwardDirection) < S(0.0))
            {
                std::reverse
                (
                    orientedFaceNodes.begin(),
                    orientedFaceNodes.end()
                );
            }

            faceConnectivity.insert
            (
                faceConnectivity.end(),
                orientedFaceNodes.begin(),
                orientedFaceNodes.end()
            );
            faceOffsets.push_back
            (
                static_cast<long long>(faceConnectivity.size())
            );
            polyhedronToFaces.push_back
            (
                static_cast<long long>(polyhedronToFaces.size())
            );
        }

        connectivity.insert
        (
            connectivity.end(),
            uniquePointIds.begin(),
            uniquePointIds.end()
        );
        offsets.push_back(static_cast<long long>(connectivity.size()));
        polyhedronOffsets.push_back
        (
            static_cast<long long>(polyhedronToFaces.size())
        );
    }

    // Only the placement is global, the values stay piece-local
    const HDF5::SlabLayout pointRows =
        HDF5::distributedRows(allNodes.size());
    const HDF5::SlabLayout cellRows = HDF5::distributedRows(numOwnedCells);
    const HDF5::SlabLayout connectivityRows =
        HDF5::distributedRows(connectivity.size());
    const HDF5::SlabLayout faceRows =
        HDF5::distributedRows(polyhedronToFaces.size());
    const HDF5::SlabLayout faceConnectivityRows =
        HDF5::distributedRows(faceConnectivity.size());

    // The NumberOf* part tables: one row per rank (length 1 when serial)
    const HDF5::SlabLayout partRow = HDF5::oneRowPerRank();

    const long long numPoints = static_cast<long long>(allNodes.size());
    const long long numCells = static_cast<long long>(numOwnedCells);
    const long long numConnectivityIds =
        static_cast<long long>(connectivity.size());
    const long long numFaces =
        static_cast<long long>(polyhedronToFaces.size());
    const long long numFaceConnectivityIds =
        static_cast<long long>(faceConnectivity.size());

    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfPoints",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numPoints,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "Points",
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        points.data(),
        pointRows,
        3
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfCells",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numCells,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfConnectivityIds",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numConnectivityIds,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "Connectivity",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        connectivity.data(),
        connectivityRows,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "Offsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        offsets.data(),
        HDF5::perPieceOffsetsRows(cellRows),
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "Types",
        H5T_STD_U8LE,
        H5T_NATIVE_UCHAR,
        types.data(),
        cellRows,
        0
    );

    // Polyhedron face datasets
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfFaces",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numFaces,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfFaceConnectivityIds",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numFaceConnectivityIds,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "NumberOfPolyhedronToFaceIds",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numFaces,
        partRow,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "FaceConnectivity",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        faceConnectivity.data(),
        faceConnectivityRows,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "FaceOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        faceOffsets.data(),
        HDF5::perPieceOffsetsRows(faceRows),
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "PolyhedronToFaces",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        polyhedronToFaces.data(),
        faceRows,
        0
    );
    HDF5::writeDataset
    (
        vtkhdfGroup_,
        "PolyhedronOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        polyhedronOffsets.data(),
        HDF5::perPieceOffsetsRows(cellRows),
        0
    );

    numCells_ = numOwnedCells;
    numNodes_ = allNodes.size();
    globalNumCells_ = static_cast<Count>(cellRows.globalRows);
    cellRowOffset_ = static_cast<Index>(cellRows.rowOffset);
    geometryWritten_ = true;

    closeHandles();

    if (debug_)
    {
        std::cout
            << "VTKHDF geometry written: " << fileName_ << '\n'
            << "  - Number of points: " << pointRows.globalRows << '\n'
            << "  - Number of cells: " << globalNumCells_ << '\n'
            << "  - Number of pieces: " << Comm::numProcessors() << '\n'
            << "  - Cell type: VTK_POLYHEDRON" << '\n';
    }
}


void HDF5CellData::appendTimeStep
(
    Scalar time,
    const ScalarFieldMap& scalarFields,
    const VectorFieldMap& vectorFields
)
{
    if (!geometryWritten_)
    {
        FatalError
        (
            "VTKHDF appendTimeStep() called before writeGeometry(): "
          + fileName_
        );
    }

    if (finalized_)
    {
        FatalError
        (
            "VTKHDF appendTimeStep() called after finalize(): " + fileName_
        );
    }

    openHandles();

    Count processedArrays = 0;

    // Scalar cell fields
    std::vector<float> scalarBuffer(numCells_);

    for (const auto& [fieldName, scalarField] : scalarFields)
    {
        if (!scalarField) continue;

        // Fields may carry a ghost tail past the owned prefix
        if (scalarField->size() < numCells_)
        {
            FatalError
            (
                "VTKHDF scalar field '" + fieldName + "' has "
              + std::to_string(scalarField->size())
              + " values but the mesh has "
              + std::to_string(numCells_) + " cells."
            );
        }

        checkFieldConsistency(fieldName, 1);

        for (Index cellIdx = 0; cellIdx < numCells_; ++cellIdx)
        {
            scalarBuffer[cellIdx] =
                static_cast<float>((*scalarField)[cellIdx]);
        }

        appendCellDataArray(fieldName, scalarBuffer.data(), 1);
        ++processedArrays;
    }

    // Vector cell fields (supplied as three scalar component fields)
    std::vector<float> vectorBuffer(3 * numCells_);

    for (const auto& [fieldName, components] : vectorFields)
    {
        if (!components[0] || !components[1] || !components[2]) continue;

        const ScalarField& cx = *components[0];
        const ScalarField& cy = *components[1];
        const ScalarField& cz = *components[2];

        if
        (
            cx.size() < numCells_
         || cy.size() < numCells_
         || cz.size() < numCells_
        )
        {
            FatalError
            (
                "VTKHDF vector field '" + fieldName
              + "' has component sizes that do not match the mesh cell "
                "count " + std::to_string(numCells_) + "."
            );
        }

        checkFieldConsistency(fieldName, 3);

        for (Index cellIdx = 0; cellIdx < numCells_; ++cellIdx)
        {
            vectorBuffer[3 * cellIdx] = static_cast<float>(cx[cellIdx]);
            vectorBuffer[3 * cellIdx + 1] = static_cast<float>(cy[cellIdx]);
            vectorBuffer[3 * cellIdx + 2] = static_cast<float>(cz[cellIdx]);
        }

        appendCellDataArray(fieldName, vectorBuffer.data(), 3);
        ++processedArrays;
    }

    if (stepCount_ > 0 && processedArrays != cellDataComponents_.size())
    {
        FatalError
        (
            "VTKHDF cell field set changed between time steps: "
          + std::to_string(processedArrays) + " arrays written but "
          + std::to_string(cellDataComponents_.size())
          + " were registered on the first step."
        );
    }

    appendStepBookkeeping(time);
    ++stepCount_;

    // Refresh NSteps, then release the file until the next write
    HDF5::writeIntScalarAttribute
    (
        stepsGroup_,
        "NSteps",
        static_cast<long long>(stepCount_)
    );
    closeHandles();
}


void HDF5CellData::finalize()
{
    if (finalized_)
    {
        return;
    }

    closeHandles();

    finalized_ = true;
}

// ****************************** Private Methods *****************************

void HDF5CellData::openHandles()
{
    file_ = HDF5::openFile(fileName_);
    vtkhdfGroup_ = HDF5::openGroup(file_, "VTKHDF");
    cellDataGroup_ = HDF5::openGroup(vtkhdfGroup_, "CellData");
    stepsGroup_ = HDF5::openGroup(vtkhdfGroup_, "Steps");
    cellDataOffsetsGroup_ = HDF5::openGroup(stepsGroup_, "CellDataOffsets");
}


void HDF5CellData::closeHandles()
{
    HDF5::closeGroup(cellDataOffsetsGroup_);
    HDF5::closeGroup(stepsGroup_);
    HDF5::closeGroup(cellDataGroup_);
    HDF5::closeGroup(vtkhdfGroup_);
    HDF5::closeFile(file_);
}


void HDF5CellData::appendCellDataArray
(
    const Name& name,
    const float* values,
    Count components
)
{
    const hsize_t columns = components == 1 ? 0 : components;

    // Chunk rows derive from the GLOBAL count for the collective create
    const HDF5::SlabLayout dataRows
    {
        static_cast<hsize_t>(numCells_),
        static_cast<hsize_t>(globalNumCells_),
        static_cast<hsize_t>(cellRowOffset_)
    };

    HDF5::appendRows
    (
        cellDataGroup_,
        name.c_str(),
        H5T_IEEE_F32LE,
        H5T_NATIVE_FLOAT,
        values,
        dataRows,
        columns,
        std::min<hsize_t>(globalNumCells_, HDF5::maxChunkRows)
    );

    const long long offsetValue =
        static_cast<long long>(cellDataRowOffsets_[name]);

    HDF5::appendRows
    (
        cellDataOffsetsGroup_,
        name.c_str(),
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &offsetValue,
        HDF5::replicatedRows(1),
        0,
        HDF5::stepChunkRows
    );

    cellDataRowOffsets_[name] += globalNumCells_;
}


void HDF5CellData::checkFieldConsistency
(
    const Name& name,
    Count components
)
{
    if (stepCount_ == 0)
    {
        cellDataComponents_[name] = components;
        cellDataRowOffsets_[name] = 0;
        return;
    }

    const auto registered = cellDataComponents_.find(name);

    if (registered == cellDataComponents_.end())
    {
        FatalError
        (
            "VTKHDF cell field '" + name
          + "' appeared after the first time step; the field set must be "
            "identical on every step."
        );
    }

    if (registered->second != components)
    {
        FatalError
        (
            "VTKHDF cell field '" + name
          + "' changed component count between time steps."
        );
    }
}


void HDF5CellData::appendStepBookkeeping(Scalar time)
{
    const double timeValue = static_cast<double>(time);
    const long long zero = 0;
    const long long numParts = static_cast<long long>(Comm::numProcessors());

    // Identical on every rank: the master writes, the others just join
    const HDF5::SlabLayout oneRow = HDF5::replicatedRows(1);

    HDF5::appendRows
    (
        stepsGroup_,
        "Values",
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        &timeValue,
        oneRow,
        0,
        HDF5::stepChunkRows
    );

    // Static geometry: every step re-reads all parts from offset zero
    HDF5::appendRows
    (
        stepsGroup_,
        "PartOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "NumberOfParts",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &numParts,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "PointOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );

    // [NSteps, 1]: one topology for an UnstructuredGrid
    HDF5::appendRows
    (
        stepsGroup_,
        "CellOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        1,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "ConnectivityIdOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        1,
        HDF5::stepChunkRows
    );

    // Polyhedron step offsets, likewise all zero
    HDF5::appendRows
    (
        stepsGroup_,
        "FaceConnectivityOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "FaceOffsetsOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
    HDF5::appendRows
    (
        stepsGroup_,
        "PolyhedronToFaceIdOffsets",
        H5T_STD_I64LE,
        H5T_NATIVE_LLONG,
        &zero,
        oneRow,
        0,
        HDF5::stepChunkRows
    );
}

} // namespace VTK
