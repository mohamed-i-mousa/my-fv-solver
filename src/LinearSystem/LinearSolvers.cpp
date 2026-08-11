/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LinearSolvers.cpp
 * @brief PETSc-backed solve workflow and runtime selection of linear solvers
 *****************************************************************************/

// ********************************** Headers *********************************

#include "LinearSolvers.h"

// Standard library headers
#include <algorithm>
#include <cmath>

// Project headers
#include "Comm.h"
#include "ErrorHandler.h"
#include "PETScRuntime.h"
#include "Reduce.h"
#include "RuntimeSelection.h"

// ************************* Special Member Functions *************************

PetscLinearSolver::PetscLinearSolver
(
    Name name,
    KSPType kspType,
    Scalar tolerance,
    Count maxIterations,
    const Name& optionsPrefix
)
:
    LinearSolver(std::move(name), tolerance, maxIterations)
{
    CheckPETSc(KSPCreate(PETScRuntime::comm(), &ksp_));
    CheckPETSc(KSPSetType(ksp_, kspType));

    // Scope case-file petscOptions entries to this equation's solver
    CheckPETSc(KSPSetOptionsPrefix(ksp_, optionsPrefix.c_str()));

    // Jacobi baseline, overridable through the case file's petscOptions
    PC preconditioner = nullptr;
    CheckPETSc(KSPGetPC(ksp_, &preconditioner));
    CheckPETSc(PCSetType(preconditioner, PCJACOBI));

    // Converge on the true residual |r|/|b|, matching the diagnostics
    CheckPETSc(KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED));

    // The current field values seed the Krylov iteration
    CheckPETSc(KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE));

    // Apply any case-file petscOptions overrides last
    CheckPETSc(KSPSetFromOptions(ksp_));
}


PetscLinearSolver::~PetscLinearSolver() noexcept
{
    if (KSPDestroy(&ksp_) != PETSC_SUCCESS)
    {
        Warning("KSPDestroy failed");
    }
    if (VecDestroy(&solution_) != PETSC_SUCCESS)
    {
        Warning("VecDestroy failed");
    }
}

// ******************************* Solver Method ******************************

void PetscLinearSolver::solve
(
    std::span<Scalar> x,
    Mat A,
    Vec b
)
{
    PetscInt systemSize = 0;
    CheckPETSc(VecGetLocalSize(b, &systemSize));

    if (static_cast<Count>(systemSize) != x.size())
    {
        FatalError("LinearSolver: x and b size mismatch");
    }

    // Solution wrapper sized on the first solve, reused afterwards
    if (solution_ == nullptr)
    {
        CheckPETSc(VecCreate(PETScRuntime::comm(), &solution_));
        CheckPETSc(VecSetSizes(solution_, systemSize, PETSC_DECIDE));
        CheckPETSc(VecSetType(solution_, VECSTANDARD));
    }

    // Keep the previous iterate for the non-finite rollback guard
    previousSolution_.assign(x.begin(), x.end());

    CheckPETSc(KSPSetOperators(ksp_, A, A));

    CheckPETSc
    (
        KSPSetTolerances
        (
            ksp_,
            tolerance(),
            PETSC_CURRENT,
            PETSC_UNLIMITED,
            static_cast<PetscInt>(maxIterations())
        )
    );

    PetscReal normB = S(0.0);
    CheckPETSc(VecNorm(b, NORM_2, &normB));

    // Zero-copy solve: the KSP writes straight into the caller's field
    CheckPETSc(VecPlaceArray(solution_, x.data()));
    CheckPETSc(KSPSolve(ksp_, b, solution_));
    CheckPETSc(VecResetArray(solution_));

    bool finite = true;
    for (const Scalar value : x)
    {
        if (!std::isfinite(value))
        {
            finite = false;
            break;
        }
    }

    // A non-finite solution on any rank must roll back every rank
    const bool anyNonFinite = globalOr(!finite);

    if (anyNonFinite)
    {
        if (Comm::master())
        {
            Warning
            (
                name()
              + ": non-finite solution, rolling back to previous iterate"
            );
        }

        std::copy(previousSolution_.begin(), previousSolution_.end(),
                  x.begin());
    }

    PetscInt iterations = 0;
    CheckPETSc(KSPGetIterationNumber(ksp_, &iterations));

    PetscReal residualNorm = S(0.0);
    CheckPETSc(KSPGetResidualNorm(ksp_, &residualNorm));

    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    CheckPETSc(KSPGetConvergedReason(ksp_, &reason));

    const SolvePerformance performance
    {
        .solverName    = name(),
        .iterations    = static_cast<int>(iterations),
        .finalResidual = normB > S(0.0) ? residualNorm / normB : residualNorm,
        .converged     = !anyNonFinite && reason > 0
    };

    setLastPerformance(performance);
}

// **************************** Runtime Selection ****************************

std::unique_ptr<LinearSolver> LinearSolver::create
(
    const Name& solverName,
    Scalar tolerance,
    Count maxIterations,
    const Name& optionsPrefix
)
{
    if (solverName == "BiCGSTAB")
    {
        return
            std::make_unique<PetscLinearSolver>
            (
                solverName,
                KSPBCGS,
                tolerance,
                maxIterations,
                optionsPrefix
            );
    }

    if (solverName == "PCG")
    {
        return
            std::make_unique<PetscLinearSolver>
            (
                solverName,
                KSPCG,
                tolerance,
                maxIterations,
                optionsPrefix
            );
    }

    RuntimeSelection::unknownSelection
    (
        "linear solver",
        solverName,
        availableSolvers()
    );
}


NameList LinearSolver::availableSolvers()
{
    return {"BiCGSTAB", "PCG"};
}
