/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file NoSlip.h
 * @brief No-slip boundary condition, an alias of FixedValue with value 0
 *
 * @class NoSlip
 * FixedValue with value 0 under the "noSlip" parse token, registered on
 * the velocity components of wall patches.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "FixedValue.h"

// ******************************* class NoSlip *******************************

class NoSlip final : public FixedValue
{
public:

// ************************* Special Member Functions *************************

    /// Construct a zero-value velocity boundary condition
    explicit NoSlip(Field field) noexcept
    :
        FixedValue{field, S(0.0)}
    {}

// ***************************** Override Methods *****************************

    /// The case-file parse token "noSlip"
    [[nodiscard]] const Name& typeName() const noexcept override;
};