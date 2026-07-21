/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryType.h
 * @brief Abstract base class for per-(patch, field) boundary conditions
 *
 * @details A BoundaryType is one boundary condition of one field on one patch.
 * It exposes its physics as per-face math, with no BC-type dispatch anywhere
 * outside this class hierarchy:
 *
 *   - faceValue() reconstructs the boundary face value phi_f, read by gradient
 *     reconstruction, Rhie-Chow, and forces.
 *   - addToDiagonal() / addToSource() return this face's owner-cell diagonal
 *     and right-hand-side contributions for matrix assembly.
 *
 * A BoundaryType carries no patch reference and no per-face state: a boundary
 * condition may be registered for a patch with no faces on this rank. The
 * manager owns per-face geometry and owner velocity and hands each method the
 * data it needs.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <iosfwd>
#include <memory>
#include <vector>

// Project headers
#include "Scalar.h"
#include "Vector.h"
#include "StringTypes.h"
#include "Field.h"

// *************************** Forward Declarations ***************************

class CaseReader;

// **************************** class BoundaryType ****************************

class BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment - Not copyable (polymorphic identity)
    BoundaryType(const BoundaryType&) = delete;
    BoundaryType& operator=(const BoundaryType&) = delete;

    /// Move constructor and assignment - Not movable (polymorphic identity)
    BoundaryType(BoundaryType&&) = delete;
    BoundaryType& operator=(BoundaryType&&) = delete;

    /// Destructor
    virtual ~BoundaryType() noexcept = default;

// **************************** Runtime Selection *****************************

    /// Create a boundary condition of the given type for the given field
    [[nodiscard]] static std::unique_ptr<BoundaryType> create
    (
        const Name& typeName,
        Field field,
        const CaseReader& patchSection
    );

    /// Case-file-selectable type names valid for the given field
    [[nodiscard]] static NameList availableTypes(Field field);

// ****************************** Public Methods ******************************

    /// The case-file parse of this type
    [[nodiscard]] virtual const Name& typeName() const noexcept = 0;

    /// The field this boundary condition applies to
    [[nodiscard]] Field field() const noexcept { return field_; }

    /// Owner-cell diagonal contribution of this boundary face
    [[nodiscard]] virtual Scalar addToDiagonal
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        const Vector& normal
    ) const = 0;

    /// Owner-cell right-hand-side contribution of this boundary face
    [[nodiscard]] virtual Scalar addToSource
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const = 0;

    /// Reconstructed boundary face value phi_f for this face
    [[nodiscard]] virtual Scalar faceValue
    (
        Scalar ownerValue,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const = 0;

// ********************************** Traits **********************************

    /// Dirichlet-like: the face value is prescribed
    [[nodiscard]] virtual bool fixesValue() const noexcept
    {
        return false;
    }

    /// The face mass flux is identically zero (Rhie-Chow flux zeroing)
    [[nodiscard]] virtual bool constrainsZeroFlux() const noexcept
    {
        return false;
    }

    /// The boundary flux receives the explicit p'-gradient correction
    [[nodiscard]] virtual bool correctsBoundaryFlux() const noexcept
    {
        return false;
    }

    /// Wall-treatment marker: the turbulence model owns the physics
    [[nodiscard]] virtual bool isWallModelled() const noexcept
    {
        return false;
    }

    /// The face value participates in the Barth-Jespersen limiter hull
    [[nodiscard]] virtual bool contributesToLimiterHull() const noexcept
    {
        return true;
    }

    /// The p' boundary condition implied by this pressure boundary condition
    [[nodiscard]] virtual std::unique_ptr<BoundaryType>
    pressureCorrectionCompanion() const;

    /// Print this boundary condition's type and parameters
    virtual void write(std::ostream& os) const = 0;

// ***************************** Protected Methods ****************************

protected:

    /// Construct with the field identity (derived classes only)
    explicit BoundaryType(Field field) noexcept
    :
        field_{field}
    {}

// ****************************** Private Members *****************************

private:

    /// The field this boundary condition applies to
    Field field_;
};
