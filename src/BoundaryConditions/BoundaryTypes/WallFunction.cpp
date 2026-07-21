/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file WallFunction.cpp
 * @brief Wall-function marker diagnostics
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "WallFunction.h"

// Standard library headers
#include <ostream>

// ***************************** Override Methods *****************************

void WallFunction::write(std::ostream& os) const
{
    os << typeName() << " (wall function)";
}
