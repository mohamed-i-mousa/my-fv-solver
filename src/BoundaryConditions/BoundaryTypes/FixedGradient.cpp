/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedGradient.cpp
 * @brief Fixed-gradient boundary coefficients and diagnostics
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "FixedGradient.h"

// Standard library headers
#include <ostream>

// ***************************** Internal Helpers *****************************

namespace
{

const Name fixedGradientToken{"fixedGradient"};

} // namespace

// ***************************** Override Methods *****************************

const Name& FixedGradient::typeName() const noexcept
{
    return fixedGradientToken;
}


Scalar FixedGradient::addToDiagonal
(
    Scalar flux,
    Scalar,
    Scalar,
    const Vector&
) const
{
    return flux;
}


Scalar FixedGradient::addToSource
(
    Scalar flux,
    Scalar GammaSf,
    Scalar,
    Scalar normalDistance,
    const Vector&,
    const Vector&
) const
{
    return GammaSf * gradient_ - (flux * gradient_ * normalDistance);
}


Scalar FixedGradient::faceValue
(
    Scalar ownerValue,
    Scalar normalDistance,
    const Vector&,
    const Vector&
) const
{
    return ownerValue + gradient_ * normalDistance;
}


std::unique_ptr<BoundaryType> FixedGradient::pressureCorrectionCompanion() const
{
    // The base pressure gradient is already imposed, so p' carries a zero
    // prescribed gradient; correctsBoundaryFlux() still feeds the explicit
    // boundary flux correction that couples pressure and velocity here
    return std::make_unique<FixedGradient>(Field::pCorr, S(0.0));
}


void FixedGradient::write(std::ostream& os) const
{
    os << typeName() << ", Gradient: " << gradient_;
}
