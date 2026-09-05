/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ConvectionScheme.cpp
 * @brief Runtime selection of convection schemes
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "ConvectionScheme.h"

// Project headers
#include "RuntimeSelection.h"
#include "Upwind.h"
#include "SecondOrderUpwind.h"
#include "CentralDifference.h"
#include "LUST.h"

// ***************************** Runtime Selection ****************************

std::unique_ptr<ConvectionScheme> ConvectionScheme::create
(
    const Name& schemeName
)
{
    if (schemeName == "Upwind")
    {
        return std::make_unique<Upwind>();
    }

    if (schemeName == "CentralDifference")
    {
        return std::make_unique<CentralDifference>();
    }

    if (schemeName == "SecondOrderUpwind")
    {
        return std::make_unique<SecondOrderUpwind>();
    }

    if (schemeName == "LUST")
    {
        return std::make_unique<LUST>();
    }

    RuntimeSelection::unknownSelection
    (
        "convection scheme",
        schemeName,
        availableSchemes()
    );
}


NameList ConvectionScheme::availableSchemes()
{
    return {"Upwind", "CentralDifference", "SecondOrderUpwind", "LUST"};
}
