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
    PETSC_CHECK(KSPCreate(PETScRuntime::comm(), &ksp_));
    PETSC_CHECK(KSPSetType(ksp_, kspType));

    // Scope case-file petscOptions entries to this equation's solver
    PETSC_CHECK(KSPSetOptionsPrefix(ksp_, optionsPrefix.c_str()));

    // Jacobi (diagonal) preconditioning: the parity baseline; runtime
    // overridable through the case file's petscOptions string
    PC preconditioner = nullptr;
    PETSC_CHECK(KSPGetPC(ksp_, &preconditioner));
    PETSC_CHECK(PCSetType(preconditioner, PCJACOBI));

    // Converge on the true residual |r|/|b|, matching the reported
    // diagnostics (the default preconditioned norm would test |P^-1 r|)
    PETSC_CHECK(KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED));

    // The current field values seed the Krylov iteration
    PETSC_CHECK(KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE));

    // Apply any case-file petscOptions overrides last
    PETSC_CHECK(KSPSetFromOptions(ksp_));
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
    PETSC_CHECK(VecGetLocalSize(b, &systemSize));

    if (static_cast<Count>(systemSize) != x.size())
    {
        FatalError("LinearSolver: x and b size mismatch");
    }

    // Lazy solution wrapper: sized on first solve (local = owned rows),
    // reused afterwards. VECSTANDARD resolves seq/mpi with the Mat; the
    // caller's owned-prefix storage is placed over it per solve
    if (solution_ == nullptr)
    {
        PETSC_CHECK(VecCreate(PETScRuntime::comm(), &solution_));
        PETSC_CHECK(VecSetSizes(solution_, systemSize, PETSC_DECIDE));
        PETSC_CHECK(VecSetType(solution_, VECSTANDARD));
    }

    // Keep the previous iterate for the non-finite rollback guard
    previousSolution_.assign(x.begin(), x.end());

    PETSC_CHECK(KSPSetOperators(ksp_, A, A));

    // Divergence abort disabled (dtol = unlimited): BiCGSTAB residuals can
    // spike transiently and recover; aborting would hand a garbage iterate
    // to the field, and the non-finite rollback guard already protects
    PETSC_CHECK
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
    PETSC_CHECK(VecNorm(b, NORM_2, &normB));

    // Zero-copy solve: the KSP writes straight into the caller's field
    PETSC_CHECK(VecPlaceArray(solution_, x.data()));
    PETSC_CHECK(KSPSolve(ksp_, b, solution_));
    PETSC_CHECK(VecResetArray(solution_));

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
    PETSC_CHECK(KSPGetIterationNumber(ksp_, &iterations));

    PetscReal residualNorm = S(0.0);
    PETSC_CHECK(KSPGetResidualNorm(ksp_, &residualNorm));

    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    PETSC_CHECK(KSPGetConvergedReason(ksp_, &reason));

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
