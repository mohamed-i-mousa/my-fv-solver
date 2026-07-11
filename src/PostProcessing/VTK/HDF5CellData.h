/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HDF5CellData.h
 * @brief VTKHDF (.vtkhdf) cell data writer for 3D polyhedral cell data
 *
 * @details Writes a static-geometry temporal VTKHDF UnstructuredGrid through
 * the HDF5 C library: the polyhedral mesh geometry is written exactly once,
 * then one time step of cell-centered fields is appended per call
 * The file is fully closed between writes, so it stays readable
 * (e.g. in ParaView) while the run continues.
 *
 * Parallel runs share ONE file: every rank holds a writer for the same
 * path and each rank's owned cells form one piece (part) of a multi-piece
 * VTKHDF layout, written collectively via MPI-IO. Every method is
 * collective — all ranks must call it, with an identical field set.
 *
 * @class HDF5CellData
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <array>
#include <map>

// External library headers
#include <hdf5.h>

// Project headers
#include "Scalar.h"
#include "Integer.h"
#include "StringTypes.h"
#include "Mesh.h"
#include "CellData.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

/// Cell-centered scalar fields keyed by output name
using ScalarFieldMap = std::map<Name, const ScalarField*>;

/// Cell-centered vector fields, supplied as three scalar component fields
using VectorFieldMap = std::map<Name, std::array<const ScalarField*, 3>>;

// **************************** class HDF5CellData ****************************

class HDF5CellData
{
public:

// ************************* Special Member Functions *************************

    /// Create the .vtkhdf file and its VTKHDF skeleton
    HDF5CellData
    (
        FilePath fileName,
        const Mesh& mesh,
        bool debug = false
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    HDF5CellData(const HDF5CellData&) = delete;
    HDF5CellData& operator=(const HDF5CellData&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    HDF5CellData(HDF5CellData&&) = delete;
    HDF5CellData& operator=(HDF5CellData&&) = delete;

    /// Destructor
    ~HDF5CellData() noexcept;

// ****************************** Public Methods ******************************

    /// Write the static polyhedral geometry once
    void writeGeometry();

    /// Append one time step of cell-centered fields
    void appendTimeStep
    (
        Scalar time,
        const ScalarFieldMap& scalarFields,
        const VectorFieldMap& vectorFields
    );

    /// Close all HDF5 handles
    void finalize();

    /// Number of time steps appended so far
    [[nodiscard]] Count stepCount() const noexcept
    {
        return stepCount_;
    }

// ****************************** Private Methods *****************************

private:

    /// Reopen the file and the VTKHDF groups for the next write
    void openHandles();

    /// Close all HDF5 handles
    void closeHandles();

    /// Append one packed cell-data buffer and its running offset row
    void appendCellDataArray
    (
        const Name& name,
        const float* values,
        Count components
    );

    /// Register (first step) or check (later steps) an array's components
    void checkFieldConsistency(const Name& name, Count components);

    /// Append the per-step time value and all-zero geometry offset rows
    void appendStepBookkeeping(Scalar time);

// ****************************** Private Members *****************************

    /// Non-owning mesh reference
    const Mesh& mesh_;

    /// Output file path
    FilePath fileName_;

    /// Verbose diagnostics flag
    bool debug_;

    /// HDF5 file handle (open only inside a write call)
    hid_t file_;

    /// /VTKHDF group
    hid_t vtkhdfGroup_;

    /// /VTKHDF/CellData group
    hid_t cellDataGroup_;

    /// /VTKHDF/Steps group
    hid_t stepsGroup_;

    /// /VTKHDF/Steps/CellDataOffsets group
    hid_t cellDataOffsetsGroup_;

    /// True once writeGeometry() has run
    bool geometryWritten_;

    /// True once finalize() has run
    bool finalized_;

    /// Number of time steps appended
    Count stepCount_;

    /// Cached owned-cell count of this rank's piece (set by writeGeometry)
    Count numCells_;

    /// Cached node count of this rank's piece (set by writeGeometry)
    Count numNodes_;

    /// Global cell count across every rank (set by writeGeometry)
    Count globalNumCells_;

    /// Global row where this rank's cell slab starts (set by writeGeometry)
    Index cellRowOffset_;

    /// Per-array cumulative cell-data row offset
    std::map<Name, Count> cellDataRowOffsets_;

    /// Per-array component count, for cross-step consistency checks
    std::map<Name, Count> cellDataComponents_;
};

} // namespace VTK
