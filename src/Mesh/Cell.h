/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Cell.h
 * @brief Represents a computational cell in the mesh
 *
 * @details This header defines the Cell class, which represents a finite
 * control volume in the computational mesh. The cell is the primary entity
 * where flow variables (pressure, velocity, etc.) are stored and solved.
 * The cell is defined by a collection of bounding faces that form a closed
 * volume.
 *
 * @class Cell
 * - Topological connectivity (bounding faces, neighboring cells)
 * - Orientation data (face alignment signs relative to cell)
 * - Cell properties (centroid, volume)
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <cstdint>
#include <vector>
#include <utility>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "Face.h"
#include "Integer.h"

// ******************************** class Cell ********************************

class Cell
{
public:

    using FaceSignList = std::vector<int8_t>;

// ************************* Special Member Functions *************************

    /// Default constructor
    Cell() noexcept = default;

    /**
     * @brief Constructs cell with connectivity data
     * @param cellIdx Unique cell identifier
     * @param faces Indices of bounding faces
     * @param neighbors Indices of neighboring cells
     * @param signs Face normal direction signs
     */
    Cell
    (
        Index cellIdx,
        IndexList faces,
        IndexList neighbors,
        FaceSignList signs
    )
    : 
        idx_(cellIdx),
        faceIndices_(std::move(faces)),
        neighborCellIndices_(std::move(neighbors)),
        faceSigns_(std::move(signs))
    {}

// ****************************** Setter Methods ******************************

    /// Set cell identifier
    void setIdx(Index cellIdx) noexcept
    {
        idx_ = cellIdx;
    }

    /// Add a bounding face with its normal direction sign
    void addFace(Index faceIdx, int8_t sign)
    {
        faceIndices_.push_back(faceIdx);
        faceSigns_.push_back(sign);
    }

    /// Set all neighbor cell indices
    void setNeighborCellIndices(const IndexList& neighbors)
    {
        neighborCellIndices_.assign(neighbors.begin(), neighbors.end());
    }

    /// Set geometry directly: ghost stubs take the owning rank's values
    void setGeometry(const Vector& centroid, Scalar volume) noexcept
    {
        centroid_ = centroid;
        volume_ = volume;
    }

// ***************************** Accessor Methods *****************************

    /// Get cell identifier
    [[nodiscard]] Index idx() const noexcept
    {
        return idx_;
    }

    /// Get bounding face indices
    [[nodiscard]] const IndexList& faceIndices() const noexcept
    {
        return faceIndices_;
    }

    /// Get neighboring cell indices
    [[nodiscard]] const IndexList& neighborCellIndices() const noexcept
    {
        return neighborCellIndices_;
    }

    /// Get face normal direction signs
    [[nodiscard]] const FaceSignList& faceSigns() const noexcept
    {
        return faceSigns_;
    }

    /// Get cell centroid
    [[nodiscard]] const Vector& centroid() const noexcept
    {
        return centroid_;
    }

    /// Get cell volume
    [[nodiscard]] Scalar volume() const noexcept
    {
        return volume_;
    }

// ************************ Geometric Property Methods ************************

    /// Calculate cell volume and centroid
    void geometricProperties
    (
        const std::vector<FaceIntegrals>& allFaceIntegrals
    );

// ****************************** Private Members *****************************

private:

    /// Unique cell identifier
    Index idx_ = 0;

    /// Indices of faces that bound this cell
    IndexList faceIndices_;

    /// Indices of neighboring cells
    IndexList neighborCellIndices_;

    /// Face normal direction signs (+1 outward, -1 inward)
    FaceSignList faceSigns_;

    /// Cell geometric center
    Vector centroid_;

    /// Cell volume
    Scalar volume_ = S(0.0);
};

// *************************** Non-Member Functions ***************************

/// Stream output operator for Cell
std::ostream& operator<<(std::ostream& os, const Cell& c);
