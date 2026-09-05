/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LUST.h
 * @brief Linear-Upwind Stabilized Transport (LUST) convection scheme
 *
 * @details LUST is a blended higher-order convection scheme combining
 * second-order central differencing with second-order linear upwind
 * differencing:
 * phi_f = alpha * phi_{f,CDS} + (1 - alpha) * phi_{f,LUD}
 * The standard blending parameter is alpha = 0.75 (75% CDS + 25% LUD). This
 * provides low numerical dissipation for scale-resolving simulations while
 * damping high-frequency dispersion errors.
 *
 * @class LUST
 * - Blended second-order convection discretization
 * - Configurable blending factor alpha (default 0.75)
 * - Deferred correction for stability within segregated momentum solves
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "ConvectionScheme.h"

// ******************************** class LUST ********************************

class LUST final : public ConvectionScheme
{
public:

// ************************* Special Member Functions *************************

    /// Constructor with blending factor (default: 0.75)
    explicit LUST(Scalar alpha = S(0.75)) noexcept
    :
        alpha_{alpha}
    {}

// ****************************** Public Methods ******************************

    /// Compute deferred correction term
    [[nodiscard]] Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;

    /// Blending factor (1 = CDS, 0 = Linear Upwind)
    [[nodiscard]] Scalar alpha() const noexcept
    {
        return alpha_;
    }

// ****************************** Private Members *****************************

private:

    /// CDS blending weight
    Scalar alpha_;
};
