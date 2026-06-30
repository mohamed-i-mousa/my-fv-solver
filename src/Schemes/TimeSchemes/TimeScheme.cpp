/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TimeScheme.cpp
 * @brief Runtime selection of time schemes
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "TimeScheme.h"

// Project headers
#include "CrankNicolson.h"
#include "ImplicitEuler.h"
#include "RuntimeSelection.h"
#include "SteadyState.h"

// **************************** Runtime Selection ****************************

std::unique_ptr<TimeScheme> TimeScheme::create
(
    const Name& schemeName,
    Scalar CrankNicolsonCoeff
)
{
    if (schemeName == "steadyState")
    {
        return std::make_unique<SteadyState>();
    }

    if (schemeName == "implicitEuler")
    {
        return std::make_unique<ImplicitEuler>();
    }

    if (schemeName == "CrankNicolson")
    {
        return std::make_unique<CrankNicolson>(CrankNicolsonCoeff);
    }

    RuntimeSelection::unknownSelection
    (
        "time scheme",
        schemeName,
        availableSchemes()
    );
}


NameList TimeScheme::availableSchemes()
{
    return {"steadyState", "implicitEuler", "CrankNicolson"};
}