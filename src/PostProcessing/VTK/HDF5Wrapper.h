/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HDF5Wrapper.h
 * @brief Thin helpers over the HDF5 C API for VTKHDF file authoring
 *
 * @details Shared by the VTKHDF volume and boundary writers: file and group
 * creation, attribute writing, fixed (write-once) datasets, and extensible
 * (append-per-time-step) datasets grown in place via hyperslab writes. Every
 * failing HDF5 status funnels into FatalError, so callers never check return
 * codes.
 *
 * @note Dataset rank is encoded by the columns parameter: 0 writes a
 * one-dimensional dataset of length rows; >= 1 writes rows x columns.
 * VTKHDF needs the distinction (Steps/Values is [NSteps] while
 * Steps/CellOffsets is [NSteps, NTopologies] even when NTopologies is 1).
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// External library headers
#include <hdf5.h>

// Project headers
#include "Integer.h"
#include "StringTypes.h"

// **************************** namespace VTK::HDF5 ***************************

namespace VTK::HDF5
{

// ******************************** Constants *********************************

/// Gzip compression level applied to chunked datasets
inline constexpr unsigned int deflateLevel = 4;

/// Chunk rows for the small per-step bookkeeping arrays (Steps group)
inline constexpr hsize_t stepChunkRows = 64;

/// Upper bound on chunk rows for large data arrays
inline constexpr hsize_t maxChunkRows = hsize_t{1} << 20;

/// Row count below which a fixed dataset stays contiguous (no compression)
inline constexpr hsize_t compressionThreshold = 1024;

// ***************************** File and Groups ******************************

/// Create (truncate) an HDF5 file and return its open handle
[[nodiscard]] hid_t createFile(const FilePath& path);

/// Open an existing HDF5 file for read-write and return its handle
[[nodiscard]] hid_t openFile(const FilePath& path);

/// Create a child group under an open location and return its handle
[[nodiscard]] hid_t createGroup(hid_t location, const char* name);

/// Open an existing group (name may be a path) and return its handle
[[nodiscard]] hid_t openGroup(hid_t location, const char* name);

/// Close a group handle and invalidate it; warns (non-fatal) on failure
void closeGroup(hid_t& group);

/// Close a file handle and invalidate it; warns (non-fatal) on failure
void closeFile(hid_t& file);

// ******************************** Attributes ********************************

/// Write a one-dimensional 64-bit integer array attribute
void writeIntArrayAttribute
(
    hid_t location,
    const char* name,
    const long long* values,
    hsize_t count
);

/// Write (or overwrite in place) a scalar 64-bit integer attribute
void writeIntScalarAttribute
(
    hid_t location,
    const char* name,
    long long value
);

/// Write a fixed-length ASCII string attribute (VTKHDF "Type" encoding)
void writeStringAttribute
(
    hid_t location,
    const char* name,
    const Name& value
);

// ********************************* Datasets *********************************

/// Write a fixed dataset in one shot (compressed above the threshold)
void writeDataset
(
    hid_t location,
    const char* name,
    hid_t fileType,
    hid_t memType,
    const void* data,
    hsize_t rows,
    hsize_t columns
);

/// Append rows to an extensible dataset
void appendRows
(
    hid_t location,
    const char* name,
    hid_t fileType,
    hid_t memType,
    const void* data,
    hsize_t rows,
    hsize_t columns,
    hsize_t chunkRows
);

} // namespace VTK::HDF5
