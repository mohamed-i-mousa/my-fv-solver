/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SecondOrderImplicit.h
 * @brief Three-time-level second-order backward implicit time scheme
 *
 * @details This header defines the SecondOrderImplicit class, which implements
 * the three-time-level second-order backward differentiation formula (BDF2).
 * For a cell of volume V and time step dt, the time derivative is discretized
 * as V/dt * (1.5 * phi^{n+1} - 2.0 * phi^n + 0.5 * phi^{n-1}). It adds
 * 1.5 * V/dt to the diagonal and V/dt * (2.0 * phi^n - 0.5 * phi^{n-1}) to the
 * source. On the first time step, it starts up using first-order implicit Euler.
 *
 * @class SecondOrderImplicit
 * - Three-time-level second-order backward implicit time integration
 * - Automatic first-order implicit Euler startup for the initial time step
 * - History tracking of phi^{n-1} through updateDdtPrevStep
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "TimeScheme.h"

// ************************ class SecondOrderImplicit *************************

class SecondOrderImplicit final : public TimeScheme
{
public:

// ************************* Special Member Functions *************************

    /// Default constructor
    SecondOrderImplicit() noexcept = default;

// ****************************** Public Methods ******************************

    /// This scheme is transient
    [[nodiscard]] bool isTransient() const noexcept override
    {
        return true;
    }

    /// Per-cell diagonal and source contributions for BDF2
    [[nodiscard]] TimeContribution coefficients
    (
        Scalar volume,
        Scalar deltaT,
        Scalar phiPrevStep,
        Scalar ddtPrevStep
    ) const noexcept override
    {
        const Scalar rDeltaT = volume / deltaT;

        if (!hasPrevPrevStep_)
        {
            // Step 1 startup: Backward Euler
            return {rDeltaT, rDeltaT * phiPrevStep};
        }

        // (1.5 * phi^{n+1} - 2.0 * phi^n + 0.5 * phi^{n-1}) / deltaT
        return
        {
            S(1.5) * rDeltaT,
            rDeltaT * (S(2.0) * phiPrevStep - S(0.5) * ddtPrevStep)
        };
    }

    /// Roll old time values: phi^n becomes the new phi^{n-1}
    [[nodiscard]] Scalar updateDdtPrevStep
    (
        Scalar /*volume*/,
        Scalar /*deltaT*/,
        Scalar /*phiNew*/,
        Scalar phiPrevStep,
        Scalar /*ddtPrevStep*/
    ) const noexcept override
    {
        hasPrevPrevStep_ = true;
        return phiPrevStep;
    }

    /// Reset startup state to initial step
    void reset() noexcept
    {
        hasPrevPrevStep_ = false;
    }

    /// Whether the scheme has advanced past the startup step
    [[nodiscard]] bool hasPrevPrevStep() const noexcept
    {
        return hasPrevPrevStep_;
    }

// ****************************** Private Members *****************************

private:

    /// Whether history field holds a valid phi^{n-1}
    mutable bool hasPrevPrevStep_ = false;
};
