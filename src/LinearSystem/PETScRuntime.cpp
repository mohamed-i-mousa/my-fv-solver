/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PETScRuntime.cpp
 * @brief RAII ownership of the PETSc runtime
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PETScRuntime.h"

// ************************* Special Member Functions *************************

PETScRuntime::PETScRuntime()
{
    CheckPETSc(PetscInitializeNoArguments());
}


PETScRuntime::~PETScRuntime() noexcept
{
    if (PetscFinalize() != PETSC_SUCCESS)
    {
        // Destructors cannot throw: report and continue shutdown
        Warning("PetscFinalize failed");
    }
}

// ****************************** Public Methods ******************************

void PETScRuntime::insertOptions(const Name& options)
{
    if (options.empty())
    {
        return;
    }

    CheckPETSc(PetscOptionsInsertString(nullptr, options.c_str()));
}
