/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedGradient.h
 * @brief Fixed-gradient (Neumann) boundary condition
 *
 * @class FixedGradient
 * Prescribes the boundary-normal gradient g: phi_f = phi_P + g * dn.
 * Coefficients a = 1, b = g * normalDistance; c = 0, d = g (the diffusive
 * flux is the specified gradient itself, Gamma_f * g * |Sf|).
 *
 * @note Case-file selectable as "fixedGradient" on any solved field (a
 * velocity component reads its entry of the "gradient" vector). On pressure
 * it drives the explicit boundary flux correction via correctsBoundaryFlux().
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// **************************** class FixedGradient ***************************

class FixedGradient final : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct with the prescribed boundary-normal gradient
    FixedGradient(Field field, Scalar gradient) noexcept
    :
        BoundaryType{field},
        gradient_{gradient}
    {}

// ***************************** Override Methods *****************************

    /// The parse token "fixedGradient"
    [[nodiscard]] const Name& typeName() const noexcept override;

    /// diag = flux (a = 1, c = 0)
    [[nodiscard]] Scalar addToDiagonal
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        const Vector& normal
    ) const override;

    /// source = GammaSf * gradient - flux * gradient * normalDistance
    [[nodiscard]] Scalar addToSource
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// phi_f = phi_P + gradient * normalDistance
    [[nodiscard]] Scalar faceValue
    (
        Scalar ownerValue,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// The boundary flux receives the explicit p'-gradient correction
    [[nodiscard]] bool correctsBoundaryFlux() const noexcept override
    {
        return true;
    }

    /// Fixed-gradient pressure implies a flux-corrected zero-gradient p'
    [[nodiscard]] std::unique_ptr<BoundaryType>
    pressureCorrectionCompanion() const override;

    /// Print type and prescribed gradient
    void write(std::ostream& os) const override;

// ****************************** Private Members *****************************

private:

    /// Prescribed boundary-normal gradient
    Scalar gradient_;
};
