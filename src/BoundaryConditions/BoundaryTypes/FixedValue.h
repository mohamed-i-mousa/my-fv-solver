/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FixedValue.h
 * @brief Fixed-value (Dirichlet) boundary condition
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// ***************************** class FixedValue *****************************

class FixedValue : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct with the prescribed boundary value
    FixedValue(Field field, Scalar value) noexcept
    :
        BoundaryType{field},
        value_{value}
    {}

// ***************************** Override Methods *****************************

    /// The case-file parse token "fixedValue"
    [[nodiscard]] const Name& typeName() const noexcept override;

    /// diag = GammaSf * diffMetric (Dirichlet: a = 0, c = -diffMetric)
    [[nodiscard]] Scalar addToDiagonal
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        const Vector& normal
    ) const override;

    /// source = value * (GammaSf * diffMetric - flux)
    [[nodiscard]] Scalar addToSource
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// The prescribed value is the face value
    [[nodiscard]] Scalar faceValue
    (
        Scalar ownerValue,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// Dirichlet: anchors the pressure level when registered on p
    [[nodiscard]] bool fixesValue() const noexcept override
    {
        return true;
    }

    /// Fixed pressure implies p' = 0
    [[nodiscard]] std::unique_ptr<BoundaryType>
    pressureCorrectionCompanion() const override;

    /// Print type and prescribed value
    void write(std::ostream& os) const override;

// ****************************** Private Members *****************************

private:

    /// Prescribed boundary value
    Scalar value_;
};
