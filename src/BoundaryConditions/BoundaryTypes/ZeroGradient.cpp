/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ZeroGradient.cpp
 * @brief Zero-gradient boundary coefficients and diagnostics
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "ZeroGradient.h"

// Standard library headers
#include <ostream>

// ***************************** Internal Helpers *****************************

namespace
{

const Name zeroGradientToken{"zeroGradient"};

} // namespace

// ***************************** Override Methods *****************************

const Name& ZeroGradient::typeName() const noexcept
{
    return zeroGradientToken;
}


Scalar ZeroGradient::addToDiagonal
(
    Scalar flux,
    Scalar,
    Scalar,
    const Vector&
) const
{
    // Zero-gradient: a = 1, c = 0
    return flux;
}


Scalar ZeroGradient::addToSource
(
    Scalar,
    Scalar,
    Scalar,
    Scalar,
    const Vector&,
    const Vector&
) const
{
    // Zero-gradient: b = 0, d = 0
    return S(0.0);
}


Scalar ZeroGradient::faceValue
(
    Scalar ownerValue,
    Scalar,
    const Vector&,
    const Vector&
) const
{
    return ownerValue;
}


std::unique_ptr<BoundaryType> ZeroGradient::pressureCorrectionCompanion() const
{
    return std::make_unique<ZeroGradient>(Field::pCorr);
}


void ZeroGradient::write(std::ostream& os) const
{
    os << typeName() << " (implies zero gradient)";
}
