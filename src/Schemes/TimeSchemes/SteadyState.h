/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SteadyState.h
 * @brief Steady-state (no time term) null-object time scheme
 *
 * @details SteadyState satisfies the TimeScheme interface for runs without a
 * transient term. It contributes nothing to the matrix, so a transport
 * equation tagged with this scheme assembles exactly as a steady equation.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "TimeScheme.h"

// **************************** class SteadyState *****************************

class SteadyState final : public TimeScheme
{
public:

// ****************************** Public Methods ******************************

    [[nodiscard]] bool isTransient() const noexcept override
    {
        return false;
    }

    [[nodiscard]] TimeContribution coefficients
    (
        Scalar /*volume*/,
        Scalar /*deltaT*/,
        Scalar /*phiOld*/,
        Scalar /*oldDdt*/
    ) const noexcept override
    {
        return {S(0.0), S(0.0)};
    }

    [[nodiscard]] Scalar updateOldDdt
    (
        Scalar /*volume*/,
        Scalar /*deltaT*/,
        Scalar /*phiNew*/,
        Scalar /*phiOld*/,
        Scalar /*oldDdt*/
    ) const noexcept override
    {
        return S(0.0);
    }
};