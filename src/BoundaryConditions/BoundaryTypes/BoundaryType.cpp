/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryType.cpp
 * @brief Shared behavior of the boundary-condition base class
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BoundaryType.h"

// Project headers
#include "ErrorHandler.h"
#include "CaseReader.h"
#include "RuntimeSelection.h"
#include "FixedGradient.h"
#include "FixedValue.h"
#include "NoSlip.h"
#include "WallFunction.h"
#include "ZeroGradient.h"

// **************************** Runtime Selection *****************************

std::unique_ptr<BoundaryType> BoundaryType::create
(
    const Name& typeName,
    Field field,
    const CaseReader& patchSection
)
{
    if (!RuntimeSelection::isKnown(typeName, availableTypes(field)))
    {
        RuntimeSelection::unknownSelection
        (
            "boundary condition type for field '"
          + Name(fieldToString(field)) + "'",
            typeName,
            availableTypes(field)
        );
    }

    if (typeName == "fixedValue")
    {
        // A velocity component picks its entry of the "value" vector
        Scalar value = S(0.0);

        switch (field)
        {
            case Field::Ux:
                value = patchSection.lookup<Vector>("value").x();
                break;
            case Field::Uy:
                value = patchSection.lookup<Vector>("value").y();
                break;
            case Field::Uz:
                value = patchSection.lookup<Vector>("value").z();
                break;
            default:
                value = patchSection.lookup<Scalar>("value");
                break;
        }

        return std::make_unique<FixedValue>(field, value);
    }

    if (typeName == "noSlip")
    {
        return std::make_unique<NoSlip>(field);
    }

    if (typeName == "zeroGradient")
    {
        return std::make_unique<ZeroGradient>(field);
    }

    if (typeName == "fixedGradient")
    {
        // A velocity component picks its entry of the "gradient" vector
        Scalar gradient = S(0.0);

        switch (field)
        {
            case Field::Ux:
                gradient = patchSection.lookup<Vector>("gradient").x();
                break;
            case Field::Uy:
                gradient = patchSection.lookup<Vector>("gradient").y();
                break;
            case Field::Uz:
                gradient = patchSection.lookup<Vector>("gradient").z();
                break;
            default:
                gradient = patchSection.lookup<Scalar>("gradient");
                break;
        }

        return std::make_unique<FixedGradient>(field, gradient);
    }

    // The remaining selectable tokens are the wall-function flavors
    return std::make_unique<WallFunction>(field, typeName);
}


NameList BoundaryType::availableTypes(Field field)
{
    switch (field)
    {
        case Field::Ux:
        case Field::Uy:
        case Field::Uz:
            return {"fixedValue", "fixedGradient", "noSlip", "zeroGradient"};

        case Field::p:
            return {"fixedValue", "fixedGradient", "zeroGradient"};

        case Field::k:
            return
                {
                    "fixedValue",
                    "fixedGradient",
                    "kWallFunction",
                    "zeroGradient"
                };

        case Field::omega:
            return
                {
                    "fixedValue",
                    "fixedGradient",
                    "omegaWallFunction",
                    "zeroGradient"
                };

        case Field::nut:
            return
                {
                    "fixedValue",
                    "fixedGradient",
                    "zeroGradient",
                    "nutWallFunction"
                };

        case Field::pCorr:
        default:
            // p' is derived from p, never case-file selectable
            return {};
    }
}

// ****************************** Public Methods ******************************

std::unique_ptr<BoundaryType> BoundaryType::pressureCorrectionCompanion() const
{
    FatalError
    (
        "Boundary condition type '" + typeName()
      + "' is registered on the pressure field but defines no pressure-"
        "correction companion. Override pressureCorrectionCompanion()."
    );
}
