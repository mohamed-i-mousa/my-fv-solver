/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TimeScheme.h
 * @brief Abstract base class for time-derivative discretization schemes
 *
 * @details This header defines the polymorphic time-integration interface for
 * the unsteady term d(phi)/dt in a transport equation. Each scheme reports
 * the per-cell contribution it makes to the implicit diagonal and to the
 * explicit right-hand side.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"

// *************************** struct TimeContribution ************************

struct TimeContribution
{
    /// Added to the cell's diagonal coefficient
    Scalar diag;

    /// Added to the cell's right-hand side
    Scalar source;
};

// ***************************** class TimeScheme *****************************

class TimeScheme
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment
    TimeScheme(const TimeScheme&) = delete;
    TimeScheme& operator=(const TimeScheme&) = delete;

    /// Move constructor and assignment
    TimeScheme(TimeScheme&&) = delete;
    TimeScheme& operator=(TimeScheme&&) = delete;

    /// Default destructor
    virtual ~TimeScheme() = default;

// **************************** Runtime Selection ****************************

    /// Construct the time scheme selected by name
    [[nodiscard]] static std::unique_ptr<TimeScheme> create
    (
        const Name& schemeName,
        Scalar CrankNicolsonCoeff = S(1.0)
    );

    /// Names of every selectable time scheme
    [[nodiscard]] static NameList availableSchemes();

// ****************************** Public Methods ******************************

    /// Whether this scheme is transient
    [[nodiscard]] virtual bool isTransient() const noexcept = 0;

    /// Per-cell implicit-diagonal and explicit-source contribution of d/dt
    [[nodiscard]] virtual TimeContribution coefficients
    (
        Scalar volume,
        Scalar deltaT,
        Scalar phiPrevStep,
        Scalar ddtPrevStep
    ) const noexcept = 0;

    /// New stored old time-derivative for the next step
    [[nodiscard]] virtual Scalar updateDdtPrevStep
    (
        Scalar volume,
        Scalar deltaT,
        Scalar phiNew,
        Scalar phiPrevStep,
        Scalar ddtPrevStep
    ) const noexcept = 0;

// ***************************** Protected Methods ****************************

protected:

    /// Default constructor
    TimeScheme() = default;
};
