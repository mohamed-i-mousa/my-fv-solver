/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Symmetry.h
 * @brief Symmetry-plane constraint boundary condition
 *
 * @class Symmetry
 * Mesh-derived constraint (PatchType::symmetry), never case-file selectable:
 * the loader stamps it on every solved field of a symmetry patch.
 *
 * On a velocity component i the value channel carries the mirror value
 * U_f = U_P - (U_P . n) n, affine in the component being solved:
 * a = 1 - n_i^2, b = -n_i * UnCross with UnCross = sum_{j != i} n_j U_j.
 * addToDiagonal() / addToSource() add the matching normal-diffusion coupling:
 * the implicit -n_i^2 * diffMetric on the diagonal and the deferred explicit
 * cross-component source -n_i * UnCross * diffMetric. On scalar fields the
 * plane is a plain zero-gradient (a = 1, b = 0, no diffusion coupling). The
 * zero-mass-flux constraint is carried by constrainsZeroFlux(), which keeps
 * the face flux - and with it every convection contribution here - zero.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// ****************************** class Symmetry ******************************

class Symmetry final : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct for the given field
    explicit Symmetry(Field field) noexcept
    :
        BoundaryType{field}
    {}

// ***************************** Override Methods *****************************

    /// The diagnostic token "symmetry" (mesh-derived, never parsed)
    [[nodiscard]] const Name& typeName() const noexcept override;

    /// Velocity: diag from normal-diffusion coupling; scalar: diag = flux
    [[nodiscard]] Scalar addToDiagonal
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        const Vector& normal
    ) const override;

    /// Velocity: source from the mirror cross-terms; scalar: source = 0
    [[nodiscard]] Scalar addToSource
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// Velocity: mirror value U_P - (U_P . n) n; scalar: owner value
    [[nodiscard]] Scalar faceValue
    (
        Scalar ownerValue,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// A symmetry plane carries zero normal mass flux
    [[nodiscard]] bool constrainsZeroFlux() const noexcept override
    {
        return true;
    }

    /// The mirror value stays out of the Barth-Jespersen hull on velocity
    /// (mirroring needs all three components; matches the legacy limiter)
    [[nodiscard]] bool contributesToLimiterHull() const noexcept override
    {
        return false;
    }

    /// Symmetry pressure implies symmetry p'
    [[nodiscard]] std::unique_ptr<BoundaryType>
    pressureCorrectionCompanion() const override;

    /// Print type
    void write(std::ostream& os) const override;
};
