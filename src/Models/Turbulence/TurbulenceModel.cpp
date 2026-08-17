/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TurbulenceModel.cpp
 * @brief Runtime selection of turbulence models
 *
 * @details Dispatches every selectable model by name, including the Laminar
 * null-object, the "no turbulence" selection. Laminar needs none of the
 * turbulence services, so its branch ignores all but the mesh and viscosity.
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "TurbulenceModel.h"

// Project headers
#include "BoundaryConditions.h"
#include "ConvectionScheme.h"
#include "GradientScheme.h"
#include "Laminar.h"
#include "LinearSolvers.h"
#include "Mesh.h"
#include "RuntimeSelection.h"
#include "kOmegaSST.h"

// **************************** Runtime Selection ****************************

std::unique_ptr<TurbulenceModel> TurbulenceModel::create
(
    const Name& modelName,
    const Mesh& mesh,
    const BoundaryConditions& bc,
    const TimeScheme& timeScheme,
    const GradientScheme& gradScheme,
    const ConvectionScheme& kScheme,
    LinearSolver& kSolver,
    const ConvectionScheme& omegaScheme,
    LinearSolver& omegaSolver,
    Scalar deltaT,
    Scalar nu,
    Scalar initialK,
    Scalar initialOmega,
    Scalar alphaK,
    Scalar alphaOmega,
    bool roughWall,
    bool debug
)
{
    if (modelName == "kOmegaSST")
    {
        return std::make_unique<kOmegaSST>
        (
            mesh,
            bc,
            timeScheme,
            gradScheme,
            kScheme,
            kSolver,
            omegaScheme,
            omegaSolver,
            deltaT,
            nu,
            initialK,
            initialOmega,
            alphaK,
            alphaOmega,
            roughWall,
            debug
        );
    }

    RuntimeSelection::unknownSelection
    (
        "turbulence model",
        Name(modelName),
        availableModels()
    );
}


NameList TurbulenceModel::availableModels()
{
    return {"Laminar", "kOmegaSST"};
}
