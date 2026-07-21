/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Symmetry.cpp
 * @brief Symmetry-plane boundary coefficients and diagnostics
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Symmetry.h"

// Standard library headers
#include <ostream>

// ***************************** Internal Helpers *****************************

namespace
{

const Name symmetryToken{"symmetry"};

} // namespace

// ***************************** Override Methods *****************************

const Name& Symmetry::typeName() const noexcept
{
    return symmetryToken;
}


Scalar Symmetry::addToDiagonal
(
    Scalar flux,
    Scalar GammaSf,
    Scalar diffMetric,
    const Vector& normal
) const
{
    const bool isVelocity =
        field() == Field::Ux || field() == Field::Uy || field() == Field::Uz;

    // Scalar fields mirror as zero-gradient: a = 1, c = 0
    if (!isVelocity)
    {
        return flux;
    }

    Scalar ni = S(0.0);

    switch (field())
    {
        case Field::Ux: ni = normal.x(); break;
        case Field::Uy: ni = normal.y(); break;
        case Field::Uz: ni = normal.z(); break;
        default: break;
    }

    // Tangential projection: a = 1 - ni^2, c = -ni^2 * diffMetric
    return GammaSf * (ni * ni * diffMetric) + flux * (S(1.0) - ni * ni);
}


Scalar Symmetry::addToSource
(
    Scalar flux,
    Scalar GammaSf,
    Scalar diffMetric,
    Scalar,
    const Vector& normal,
    const Vector& ownerVelocity
) const
{
    const bool isVelocity =
        field() == Field::Ux || field() == Field::Uy || field() == Field::Uz;

    // Scalar fields carry no tangential source: b = 0, d = 0
    if (!isVelocity)
    {
        return S(0.0);
    }

    Scalar ni = S(0.0);
    Scalar UnCross = S(0.0);

    switch (field())
    {
        case Field::Ux:
            ni = normal.x();
            UnCross = normal.y() * ownerVelocity.y()
                    + normal.z() * ownerVelocity.z();
            break;
        case Field::Uy:
            ni = normal.y();
            UnCross = normal.x() * ownerVelocity.x()
                    + normal.z() * ownerVelocity.z();
            break;
        case Field::Uz:
            ni = normal.z();
            UnCross = normal.x() * ownerVelocity.x()
                    + normal.y() * ownerVelocity.y();
            break;
        default:
            break;
    }

    // b = -ni * UnCross, d = b * diffMetric
    return GammaSf * (-ni * UnCross * diffMetric) - flux * (-ni * UnCross);
}


Scalar Symmetry::faceValue
(
    Scalar ownerValue,
    Scalar,
    const Vector& normal,
    const Vector& ownerVelocity
) const
{
    const bool isVelocity =
        field() == Field::Ux || field() == Field::Uy || field() == Field::Uz;

    // Scalar fields see the plane as zero gradient
    if (!isVelocity)
    {
        return ownerValue;
    }

    Scalar ni = S(0.0);
    Scalar UnCross = S(0.0);

    switch (field())
    {
        case Field::Ux:
            ni = normal.x();
            UnCross = normal.y() * ownerVelocity.y()
                    + normal.z() * ownerVelocity.z();
            break;
        case Field::Uy:
            ni = normal.y();
            UnCross = normal.x() * ownerVelocity.x()
                    + normal.z() * ownerVelocity.z();
            break;
        case Field::Uz:
            ni = normal.z();
            UnCross = normal.x() * ownerVelocity.x()
                    + normal.y() * ownerVelocity.y();
            break;
        default:
            break;
    }

    // Mirror value: U_f = (1 - ni^2) * U_P,i - ni * UnCross
    return (S(1.0) - ni * ni) * ownerValue + (-ni * UnCross);
}


std::unique_ptr<BoundaryType> Symmetry::pressureCorrectionCompanion() const
{
    return std::make_unique<Symmetry>(Field::pCorr);
}


void Symmetry::write(std::ostream& os) const
{
    os << typeName() << " (symmetry plane)";
}
