/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ZeroGradient.h
 * @brief Zero-normal-gradient boundary condition
 *
 * @class ZeroGradient
 * The face value is the owner-cell value and the diffusive flux vanishes:
 * a = 1, b = 0; c = 0, d = 0.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "BoundaryType.h"

// **************************** class ZeroGradient ****************************

class ZeroGradient : public BoundaryType
{
public:

// ************************* Special Member Functions *************************

    /// Construct for the given field
    explicit ZeroGradient(Field field) noexcept
    :
        BoundaryType{field}
    {}

// ***************************** Override Methods *****************************

    /// The case-file parse token "zeroGradient"
    [[nodiscard]] const Name& typeName() const noexcept override;

    /// diag = flux (a = 1, c = 0)
    [[nodiscard]] Scalar addToDiagonal
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        const Vector& normal
    ) const override;

    /// source = 0 (b = 0, d = 0)
    [[nodiscard]] Scalar addToSource
    (
        Scalar flux,
        Scalar GammaSf,
        Scalar diffMetric,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// The owner-cell value is the face value
    [[nodiscard]] Scalar faceValue
    (
        Scalar ownerValue,
        Scalar normalDistance,
        const Vector& normal,
        const Vector& ownerVelocity
    ) const override;

    /// Zero-gradient pressure implies zero-gradient p'
    [[nodiscard]] std::unique_ptr<BoundaryType>
    pressureCorrectionCompanion() const override;

    /// Print type
    void write(std::ostream& os) const override;
};
