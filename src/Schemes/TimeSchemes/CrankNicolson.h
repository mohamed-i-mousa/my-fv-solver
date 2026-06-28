/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CrankNicolson.h
 * @brief Second-order Crank-Nicolson time scheme
 *
 * @details Crank-Nicolson is implemented in the form that keeps
 * the spatial operator fully implicit and carries a per-cell old time
 * derivative (ddt0) between steps. With coefft = 1 + oc (oc in [0, 1]):
 *   diag   += coefft * V/dt
 *   source += coefft * V/dt * phi^n + oc * ddt0
 *   ddt0_new = coefft * V/dt * (phi^{n+1} - phi^n) - oc * ddt0
 * Because the stored ddt0 equals the previous spatial residual, oc = 1 yields
 * true second-order Crank-Nicolson and oc = 0 degenerates to backward Euler.
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
    explicit CrankNicolson(Scalar offCenteringCoefficient) noexcept
    :
        offCenteringCoeff_{offCenteringCoefficient}
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
        Scalar phiOld,
        Scalar oldDdt
    ) const noexcept override
    {
        const Scalar coefft = S(1.0) + offCenteringCoeff_;
        const Scalar rDeltaT = coefft * volume / deltaT;
        return {rDeltaT, rDeltaT * phiOld + offCenteringCoeff_ * oldDdt};
    }

    [[nodiscard]] Scalar updateOldDdt
    (
        Scalar volume,
        Scalar deltaT,
        Scalar phiNew,
        Scalar phiOld,
        Scalar oldDdt
    ) const noexcept override
    {
        const Scalar coefft = S(1.0) + offCenteringCoeff_;
        return
            coefft * (volume / deltaT) * (phiNew - phiOld)
          - offCenteringCoeff_ * oldDdt;
    }

// ****************************** Private Members *****************************

private:

    /// Off-centering coefficient (0 = backward Euler, 1 = Crank-Nicolson)
    Scalar offCenteringCoeff_;
};