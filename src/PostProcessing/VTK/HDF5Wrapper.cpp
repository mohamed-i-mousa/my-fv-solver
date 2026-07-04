/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HDF5Wrapper.cpp
 * @brief Implementation of the HDF5 C API helpers for VTKHDF authoring
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "HDF5Wrapper.h"

// Standard library headers
#include <algorithm>

// Project headers
#include "ErrorHandler.h"

// **************************** namespace VTK::HDF5 ***************************

namespace VTK::HDF5
{

// ***************************** Internal Helpers *****************************

namespace
{

// RAII over a short-lived HDF5 identifier (dataspace, plist, dataset, ...)
class ScopedHandle
{
public:

    ScopedHandle(hid_t handle, herr_t (*closer)(hid_t))
    :
        handle_{handle},
        closer_{closer}
    {}

    // Copy constructor and assignment — Not copyable (owning handle)
    ScopedHandle(const ScopedHandle&) = delete;
    ScopedHandle& operator=(const ScopedHandle&) = delete;

    // Move constructor and assignment — Not movable (owning handle)
    ScopedHandle(ScopedHandle&&) = delete;
    ScopedHandle& operator=(ScopedHandle&&) = delete;

    // Destructor — closes the handle if valid
    ~ScopedHandle() noexcept
    {
        if (handle_ >= 0)
        {
            closer_(handle_);
        }
    }

    [[nodiscard]] hid_t get() const noexcept
    {
        return handle_;
    }

private:

    hid_t handle_;
    herr_t (*closer_)(hid_t);
};


// Abort with context when an HDF5 call returns an invalid identifier
hid_t requireValid(hid_t handle, const Message& what)
{
    if (handle < 0)
    {
        FatalError("HDF5 error: failed to " + what + ".");
    }

    return handle;
}


// Abort with context when an HDF5 call returns a failing status
void requireStatus(herr_t status, const Message& what)
{
    if (status < 0)
    {
        FatalError("HDF5 error: failed to " + what + ".");
    }
}

} // namespace

// ***************************** File and Groups ******************************

hid_t createFile(const FilePath& path)
{
    return requireValid
    (
        H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
        "create file '" + path + "'"
    );
}


hid_t openFile(const FilePath& path)
{
    return requireValid
    (
        H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT),
        "open file '" + path + "'"
    );
}


hid_t createGroup(hid_t location, const char* name)
{
    return requireValid
    (
        H5Gcreate2(location, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        "create group '" + Message{name} + "'"
    );
}


hid_t openGroup(hid_t location, const char* name)
{
    return requireValid
    (
        H5Gopen2(location, name, H5P_DEFAULT),
        "open group '" + Message{name} + "'"
    );
}

void closeGroup(hid_t& group)
{
    if (group >= 0)
    {
        if (H5Gclose(group) < 0)
        {
            Warning("Failed to close an HDF5 group handle.");
        }

        group = H5I_INVALID_HID;
    }
}


void closeFile(hid_t& file)
{
    if (file >= 0)
    {
        if (H5Fclose(file) < 0)
        {
            Warning("Failed to close an HDF5 file handle.");
        }

        file = H5I_INVALID_HID;
    }
}

// ******************************** Attributes ********************************

void writeIntArrayAttribute
(
    hid_t location,
    const char* name,
    const long long* values,
    hsize_t count
)
{
    const Message what = "write attribute '" + Message{name} + "'";

    const ScopedHandle space
    {
        requireValid(H5Screate_simple(1, &count, nullptr), what),
        H5Sclose
    };

    const ScopedHandle attribute
    {
        requireValid
        (
            H5Acreate2
            (
                location,
                name,
                H5T_STD_I64LE,
                space.get(),
                H5P_DEFAULT,
                H5P_DEFAULT
            ),
            what
        ),
        H5Aclose
    };

    requireStatus(H5Awrite(attribute.get(), H5T_NATIVE_LLONG, values), what);
}


void writeIntScalarAttribute
(
    hid_t location,
    const char* name,
    long long value
)
{
    const Message what = "write attribute '" + Message{name} + "'";

    const htri_t exists = H5Aexists(location, name);
    requireStatus(exists < 0 ? -1 : 0, what);

    if (exists > 0)
    {
        const ScopedHandle attribute
        {
            requireValid(H5Aopen(location, name, H5P_DEFAULT), what),
            H5Aclose
        };

        requireStatus
        (
            H5Awrite(attribute.get(), H5T_NATIVE_LLONG, &value),
            what
        );
        return;
    }

    const ScopedHandle space
    {
        requireValid(H5Screate(H5S_SCALAR), what),
        H5Sclose
    };

    const ScopedHandle attribute
    {
        requireValid
        (
            H5Acreate2
            (
                location,
                name,
                H5T_STD_I64LE,
                space.get(),
                H5P_DEFAULT,
                H5P_DEFAULT
            ),
            what
        ),
        H5Aclose
    };

    requireStatus(H5Awrite(attribute.get(), H5T_NATIVE_LLONG, &value), what);
}


void writeStringAttribute
(
    hid_t location,
    const char* name,
    const Name& value
)
{
    const Message what = "write attribute '" + Message{name} + "'";

    // Fixed-length ASCII, the encoding VTKHDF readers expect for "Type"
    const ScopedHandle stringType
    {
        requireValid(H5Tcopy(H5T_C_S1), what),
        H5Tclose
    };

    requireStatus(H5Tset_size(stringType.get(), value.size()), what);
    requireStatus(H5Tset_strpad(stringType.get(), H5T_STR_NULLPAD), what);
    requireStatus(H5Tset_cset(stringType.get(), H5T_CSET_ASCII), what);

    const ScopedHandle space
    {
        requireValid(H5Screate(H5S_SCALAR), what),
        H5Sclose
    };

    const ScopedHandle attribute
    {
        requireValid
        (
            H5Acreate2
            (
                location,
                name,
                stringType.get(),
                space.get(),
                H5P_DEFAULT,
                H5P_DEFAULT
            ),
            what
        ),
        H5Aclose
    };

    requireStatus
    (
        H5Awrite(attribute.get(), stringType.get(), value.data()),
        what
    );
}

// ********************************* Datasets *********************************

void writeDataset
(
    hid_t location,
    const char* name,
    hid_t fileType,
    hid_t memType,
    const void* data,
    hsize_t rows,
    hsize_t columns
)
{
    const Message what = "write dataset '" + Message{name} + "'";

    const int rank = columns == 0 ? 1 : 2;
    const hsize_t dims[2] = {rows, columns};

    const ScopedHandle space
    {
        requireValid(H5Screate_simple(rank, dims, nullptr), what),
        H5Sclose
    };

    const ScopedHandle plist
    {
        requireValid(H5Pcreate(H5P_DATASET_CREATE), what),
        H5Pclose
    };

    if (rows >= compressionThreshold)
    {
        const hsize_t chunk[2] =
        {
            std::min(rows, maxChunkRows),
            columns == 0 ? 1 : columns
        };

        requireStatus(H5Pset_chunk(plist.get(), rank, chunk), what);
        requireStatus(H5Pset_shuffle(plist.get()), what);
        requireStatus(H5Pset_deflate(plist.get(), deflateLevel), what);
    }

    const ScopedHandle dataset
    {
        requireValid
        (
            H5Dcreate2
            (
                location,
                name,
                fileType,
                space.get(),
                H5P_DEFAULT,
                plist.get(),
                H5P_DEFAULT
            ),
            what
        ),
        H5Dclose
    };

    if (rows > 0)
    {
        requireStatus
        (
            H5Dwrite(dataset.get(), memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
            what
        );
    }
}


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
)
{
    const Message what = "append to dataset '" + Message{name} + "'";

    const int rank = columns == 0 ? 1 : 2;

    const htri_t exists = H5Lexists(location, name, H5P_DEFAULT);
    requireStatus(exists < 0 ? -1 : 0, what);

    if (exists == 0)
    {
        const hsize_t initial[2] = {0, columns};
        const hsize_t maximum[2] = {H5S_UNLIMITED, columns};
        const hsize_t chunk[2] =
        {
            std::max<hsize_t>(chunkRows, 1),
            columns == 0 ? 1 : columns
        };

        const ScopedHandle space
        {
            requireValid(H5Screate_simple(rank, initial, maximum), what),
            H5Sclose
        };

        const ScopedHandle plist
        {
            requireValid(H5Pcreate(H5P_DATASET_CREATE), what),
            H5Pclose
        };

        requireStatus(H5Pset_chunk(plist.get(), rank, chunk), what);
        requireStatus(H5Pset_shuffle(plist.get()), what);
        requireStatus(H5Pset_deflate(plist.get(), deflateLevel), what);

        const ScopedHandle created
        {
            requireValid
            (
                H5Dcreate2
                (
                    location,
                    name,
                    fileType,
                    space.get(),
                    H5P_DEFAULT,
                    plist.get(),
                    H5P_DEFAULT
                ),
                what
            ),
            H5Dclose
        };
    }

    if (rows == 0)
    {
        return;
    }

    const ScopedHandle dataset
    {
        requireValid(H5Dopen2(location, name, H5P_DEFAULT), what),
        H5Dclose
    };

    hsize_t current[2] = {0, 0};
    {
        const ScopedHandle currentSpace
        {
            requireValid(H5Dget_space(dataset.get()), what),
            H5Sclose
        };

        requireStatus
        (
            H5Sget_simple_extent_dims(currentSpace.get(), current, nullptr),
            what
        );
    }

    const hsize_t extended[2] = {current[0] + rows, columns};
    requireStatus(H5Dset_extent(dataset.get(), extended), what);

    const ScopedHandle fileSpace
    {
        requireValid(H5Dget_space(dataset.get()), what),
        H5Sclose
    };

    const hsize_t offset[2] = {current[0], 0};
    const hsize_t count[2] = {rows, columns};

    requireStatus
    (
        H5Sselect_hyperslab
        (
            fileSpace.get(),
            H5S_SELECT_SET,
            offset,
            nullptr,
            count,
            nullptr
        ),
        what
    );

    const ScopedHandle memSpace
    {
        requireValid(H5Screate_simple(rank, count, nullptr), what),
        H5Sclose
    };

    requireStatus
    (
        H5Dwrite
        (
            dataset.get(),
            memType,
            memSpace.get(),
            fileSpace.get(),
            H5P_DEFAULT,
            data
        ),
        what
    );
}

} // namespace VTK::HDF5
