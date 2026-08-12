/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CrankNicolson.h
 * @brief Second-order Crank-Nicolson time scheme
 *
 * @details Crank-Nicolson is implemented in a form that keeps the spatial
 * term fully implicit by carrying a per-cell previous time step derivative
 * (ddtPrevStep) between steps.
 *
 * For blending coefficient c = CrankNicolsonCoeff in [0, 1]
 * (where theta = 1 / (1 + c)):
 *   (1 + c) * V/dt * (phi^{n+1} - phi^n) = ddt^{n+1} + c * ddt^n
 *
 * Discretization terms:
 *   diag   += (1 + c) * V/dt
 *   source += (1 + c) * V/dt * phi^n + c * ddtPrevStep
 *   ddtPrevStep_new = (1 + c) * V/dt * (phi^{n+1} - phi^n) - c * ddtPrevStep
 *
 * Key Special Cases (CrankNicolsonCoeff):
 *   - c = 1.0 (theta = 0.50): Pure 2nd-order Crank-Nicolson scheme.
 *   - c = 0.9 (theta ≈ 0.53): Off-centered Crank-Nicolson
 *                             (damps high-frequency oscillations).
 *   - c = 0.0 (theta = 1.00): Degenerates to 1st-order Implicit Euler.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "TimeScheme.h"

// **************************** class CrankNicolson ***************************

class CrankNicolson final : public TimeScheme
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    explicit CrankNicolson(Scalar CrankNicolsonCoeff) noexcept
    :
        CrankNicolsonCoeff_{CrankNicolsonCoeff}
    {}

// ****************************** Public Methods ******************************

    [[nodiscard]] bool isTransient() const noexcept override
    {
        return true;
    }

    [[nodiscard]] TimeContribution coefficients
    (
        Scalar volume,
        Scalar deltaT,
        Scalar phiPrevStep,
        Scalar ddtPrevStep
    ) const noexcept override
    {
        const Scalar coefft = S(1.0) + CrankNicolsonCoeff_;
        const Scalar rDeltaT = coefft * volume / deltaT;
        
        return
            {
                rDeltaT,
                rDeltaT * phiPrevStep + CrankNicolsonCoeff_ * ddtPrevStep
            };
    }
    
    [[nodiscard]] Scalar updateDdtPrevStep
    (
        Scalar volume,
        Scalar deltaT,
        Scalar phiNew,
        Scalar phiPrevStep,
        Scalar ddtPrevStep
    ) const noexcept override
    {
        const Scalar coefft = S(1.0) + CrankNicolsonCoeff_;
        return
            coefft * (volume / deltaT) * (phiNew - phiPrevStep)
          - CrankNicolsonCoeff_ * ddtPrevStep;
    }

// ****************************** Private Members *****************************

private:

    /// Crank-Nicolson coefficient (0 = backward Euler, 1 = Crank-Nicolson)
    Scalar CrankNicolsonCoeff_;
};