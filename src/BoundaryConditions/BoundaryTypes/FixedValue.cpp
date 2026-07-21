/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedValue.cpp
 * @brief Fixed-value boundary coefficients and diagnostics
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "FixedValue.h"

// Standard library headers
#include <ostream>

// ***************************** Internal Helpers *****************************

namespace
{

const Name fixedValueToken{"fixedValue"};

} // namespace

// ***************************** Override Methods *****************************

const Name& FixedValue::typeName() const noexcept
{
    return fixedValueToken;
}


Scalar FixedValue::addToDiagonal
(
    Scalar,
    Scalar GammaSf,
    Scalar diffMetric,
    const Vector&
) const
{
    // Dirichlet: a = 0, c = -diffMetric
    return GammaSf * diffMetric;
}


Scalar FixedValue::addToSource
(
    Scalar flux,
    Scalar GammaSf,
    Scalar diffMetric,
    Scalar,
    const Vector&,
    const Vector&
) const
{
    // Dirichlet: b = value, d = value * diffMetric
    return GammaSf * (value_ * diffMetric) - flux * value_;
}


Scalar FixedValue::faceValue
(
    Scalar,
    Scalar,
    const Vector&,
    const Vector&
) const
{
    return value_;
}


std::unique_ptr<BoundaryType> FixedValue::pressureCorrectionCompanion() const
{
    return std::make_unique<FixedValue>(Field::pCorr, S(0.0));
}


void FixedValue::write(std::ostream& os) const
{
    os << typeName() << ", Value: " << value_;
}
