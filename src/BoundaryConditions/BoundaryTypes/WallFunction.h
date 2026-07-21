/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file WallFunction.h
 * @brief Wall-function marker boundary condition for turbulence fields
 *
 * @class WallFunction
 * One class covers the kWallFunction/omegaWallFunction/nutWallFunction
 * flavors, distinguished by the stored parse token. Its assembly numerics
 * are inherited zero-gradient; the wall physics (omega wall values, k
 * production override, nut log-law) lives in the turbulence model, which
 * activates it wherever isWallModelled() is true. A flavor graduates to its
 * own class only when it gains distinct boundary-layer behavior.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <utility>

// Project headers
#include "ZeroGradient.h"

// **************************** class WallFunction ****************************

class WallFunction final : public ZeroGradient
{
public:

// ************************* Special Member Functions *************************

    /// Construct with the flavor's parse token (e.g. "kWallFunction")
    WallFunction(Field field, Name typeName)
    :
        ZeroGradient{field},
        typeName_{std::move(typeName)}
    {}

// ***************************** Override Methods *****************************

    /// The flavor's case-file parse token
    [[nodiscard]] const Name& typeName() const noexcept override
    {
        return typeName_;
    }

    /// Wall-treatment marker: the turbulence model owns the physics
    [[nodiscard]] bool isWallModelled() const noexcept override
    {
        return true;
    }

    /// Print type
    void write(std::ostream& os) const override;

// ****************************** Private Members *****************************

private:

    /// The flavor's case-file parse token
    Name typeName_;
};
