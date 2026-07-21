/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file NoSlip.cpp
 * @brief No-slip boundary condition implementation
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "NoSlip.h"

// ***************************** Internal Helpers *****************************

namespace
{

const Name noSlipToken{"noSlip"};

} // namespace

// ***************************** Override Methods *****************************

const Name& NoSlip::typeName() const noexcept
{
    return noSlipToken;
}