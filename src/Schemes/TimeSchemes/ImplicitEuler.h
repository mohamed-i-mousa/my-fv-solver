/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ImplicitEuler.h
 * @brief First-order implicit (backward) Euler time scheme
 *
 * @details Backward Euler discretizes d(phi)/dt over a cell of volume V as
 * V/dt * (phi^{n+1} - phi^n). It adds V/dt to the diagonal and V/dt * phi^n
 * to the source. It is unconditionally stable and carries no extra stored
 * state.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "TimeScheme.h"

// **************************** class ImplicitEuler ***************************

class ImplicitEuler final : public TimeScheme
{
public:

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
        Scalar /*ddtPrevStep*/
    ) const noexcept override
    {
        const Scalar rDeltaT = volume / deltaT;
        return {rDeltaT, rDeltaT * phiPrevStep};
    }

    [[nodiscard]] Scalar updateDdtPrevStep
    (
        Scalar /*volume*/,
        Scalar /*deltaT*/,
        Scalar /*phiNew*/,
        Scalar /*phiPrevStep*/,
        Scalar /*ddtPrevStep*/
    ) const noexcept override
    {
        return S(0.0);
    }
};