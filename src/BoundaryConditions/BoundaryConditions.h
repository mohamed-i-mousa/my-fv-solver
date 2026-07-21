/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryConditions.h
 * @brief Manages and owns boundary conditions for the CFD solver
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <array>
#include <map>
#include <memory>
#include <vector>

// Project headers
#include "Scalar.h"
#include "MeshContainers.h"
#include "Face.h"
#include "BoundaryPatch.h"
#include "CellData.h"
#include "Field.h"
#include "Integer.h"
#include "StringTypes.h"
#include "BoundaryType.h"

// ************************* class BoundaryConditions *************************

class BoundaryConditions
{
public:

    using BoundaryTypeMap =
        std::map<Name, std::map<Field, std::unique_ptr<BoundaryType>>>;

// ****************************** Setter Methods ******************************

    /// Add a boundary patch from mesh reader
    void addPatch(BoundaryPatch patch);

    /// Register a boundary condition object for a field on a patch
    void setBoundaryType
    (
        const Name& patchName,
        Field field,
        std::unique_ptr<BoundaryType> bc
    );

    /// Seal registration and build the per-face trait flag arrays
    void finalize();

    /// Snapshot owner-cell velocity per boundary face so single-component
    /// consumers can reconstruct the symmetry mirror (rank-local)
    void snapshotBoundaryVelocity
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    );

// ***************************** Accessor Methods *****************************

    /// Get boundary patch by name
    [[nodiscard]] const BoundaryPatch& patch(const Name& patchName) const;

    /// Get all boundary patches
    [[nodiscard]] const PatchList& patches() const noexcept
    {
        return patches_;
    }

    /// Get number of patches
    [[nodiscard]] Count numPatches() const noexcept
    {
        return patches_.size();
    }

    /// Compact boundary index of a boundary face (FatalError on interior)
    [[nodiscard]] Index boundaryIdx(Index faceIdx) const;

    /// Whether a boundary condition is registered for the field at the face
    [[nodiscard]] bool isRegistered
    (
        Field field,
        Index boundaryIdx
    ) const noexcept
    {
        return boundaryTypeAt_[fieldSlot(field)][boundaryIdx] != nullptr;
    }

    /// The boundary condition object for a field at a compact boundary face
    [[nodiscard]] const BoundaryType& boundaryType
    (
        Field field,
        Index boundaryIdx
    ) const;

    /// Over-relaxed orthogonal diffusion metric of a boundary face
    [[nodiscard]] Scalar diffMetric(Index boundaryIdx) const noexcept
    {
        return diffMetric_[boundaryIdx];
    }

    /// Owner-to-face normal distance of a boundary face
    [[nodiscard]] Scalar normalDistance(Index boundaryIdx) const noexcept
    {
        return normalDistance_[boundaryIdx];
    }

    /// Unit outward normal of a boundary face
    [[nodiscard]] const Vector& normal(Index boundaryIdx) const noexcept
    {
        return normals_[boundaryIdx];
    }

    /// Owner-cell velocity snapshot of a boundary face
    [[nodiscard]] const Vector& ownerVelocity(Index boundaryIdx) const noexcept
    {
        return ownerVelocity_[boundaryIdx];
    }

    /// Whether the face's value stays out of the velocity limiter hull
    [[nodiscard]] bool excludedFromVelocityHull(Index faceIdx) const noexcept
    {
        return velocityHullExcluded_[faceIdx] != 0;
    }

    /// Whether the face carries an identically zero mass flux
    [[nodiscard]] bool constrainsZeroFlux(Index faceIdx) const noexcept
    {
        return fluxConstrained_[faceIdx] != 0;
    }

    /// Whether the face flux receives the explicit p'-gradient correction
    [[nodiscard]] bool correctsBoundaryFlux(Index faceIdx) const noexcept
    {
        return correctsFlux_[faceIdx] != 0;
    }

    /// Whether a boundary condition object is registered for the field/patch
    [[nodiscard]] bool hasBoundaryType
    (
        const Name& patchName,
        Field field
    ) const noexcept;

    /// Get the boundary condition object for a field on a patch
    [[nodiscard]] const BoundaryType& boundaryType
    (
        const Name& patchName,
        Field field
    ) const;

    /// Link boundary faces to their owning patches
    void linkFaces(FaceList& faces);

    /// Validate boundary condition patch names against mesh patch names
    void validatePatchNames() const;

    /// Print summary of all boundary conditions
    void printSummary() const;

// ****************************** Private Methods *****************************

private:

    /// Coefficient-array slot of a field
    [[nodiscard]] static Count fieldSlot(Field field) noexcept
    {
        return static_cast<Count>(field);
    }

// ****************************** Private Members *****************************

private:

    /// Number of solver fields with boundary-coefficient storage
    static constexpr Count numFields_ = static_cast<Count>(Field::nut) + 1;

    /// Sentinel marking a face or patch outside the compact boundary range
    static constexpr Index noBoundaryIdx_ = static_cast<Index>(-1);

    /// Nested map: patch name → field → boundary condition object
    BoundaryTypeMap boundaryTypes_;

    /// Per-field boundary type per compact boundary face (nullptr if none)
    std::array<std::vector<const BoundaryType*>, numFields_> boundaryTypeAt_;

    /// Shared boundary geometry cache (compact-indexed)
    IndexList geomOwnerCells_;
    std::vector<Vector> normals_;
    ScalarList diffMetric_;
    ScalarList normalDistance_;

    /// Owner-cell velocity snapshot per compact boundary face (symmetry)
    std::vector<Vector> ownerVelocity_;

    /// Global face index → compact boundary index (sentinel for interior)
    IndexList boundaryIdx_;

    /// Compact slice start per patch (sentinel for processor patches)
    IndexList patchStart_;

    /// Per-face trait flags (global face indexed, built by finalize())
    std::vector<char> fluxConstrained_;
    std::vector<char> correctsFlux_;
    std::vector<char> velocityHullExcluded_;

    /// All boundary patches
    PatchList patches_;

    /// True after linkFaces() to prevent addPatch() after linking
    bool linked_ = false;

    /// True after finalize() to seal registration
    bool finalized_ = false;
};
