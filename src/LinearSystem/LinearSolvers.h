/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LinearSolvers.h
 * @brief Iterative solver hierarchy for sparse linear systems
 *
 * @details Polymorphic LinearSolver interface backed by PETSc's KSP Krylov
 * solvers. Each PetscLinearSolver owns one KSP configured at construction
 * (solver type, Jacobi preconditioner, unpreconditioned residual norm);
 * create() maps each selectable name to its KSP type — "BiCGSTAB" to
 * KSPBCGS for the non-symmetric momentum/turbulence systems, "PCG" to
 * KSPCG for the symmetric positive definite p' system. The solution is
 * bound zero-copy over the caller's field storage, so the solve writes the
 * result straight into the field.
 *
 * Defaults set here can be overridden at runtime through the case file's
 * linearSolvers.petscOptions string (see PETScRuntime::insertOptions).
 * Each KSP reads the options database through its own equation prefix
 * (momentum_, pressure_, k_, omega_), so an option targets one solver —
 * e.g. "-pressure_pc_type icc" swaps only the pressure preconditioner
 * without a rebuild; an unprefixed "-pc_type" matches no solver.
 *
 * @class LinearSolver
 * - Abstract base with virtual solve()
 * - Holds shared convergence configuration and cached solve diagnostics
 *
 * @class PetscLinearSolver
 * - Owns the KSP and the reusable solution Vec wrapper
 * - Implements the common solve workflow once
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <limits>
#include <memory>
#include <span>
#include <utility>

// External library headers
#include <petscksp.h>

// Project headers
#include "Integer.h"
#include "Scalar.h"
#include "StringTypes.h"

// ************************** struct SolvePerformance *************************

struct SolvePerformance
{
    /// Solver name
    Name solverName = {};

    /// Iterations performed by the solve call
    int iterations = 0;

    /// Final relative residual |r|/|b| reported by the solve
    Scalar finalResidual = S(0.0);

    /// Whether the solver reported successful convergence
    bool converged = false;
};

// **************************** class LinearSolver ****************************

class LinearSolver
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    LinearSolver
    (
        Name name,
        Scalar tolerance = S(1e-6),
        Count maxIterations = 1000
    )
    :
        name_{std::move(name)},
        tolerance_{tolerance},
        maxIterations_{maxIterations}
    {}


    /// Copy constructor and assignment - Not copyable (Polymorphic)
    LinearSolver(const LinearSolver&) = delete;
    LinearSolver& operator=(const LinearSolver&) = delete;

    /// Move constructor and assignment - Not movable
    LinearSolver(LinearSolver&&) = delete;
    LinearSolver& operator=(LinearSolver&&) = delete;

    /// Virtual destructor for polymorphic deletion
    virtual ~LinearSolver() noexcept = default;

// **************************** Runtime Selection ****************************

    /// Construct the solver by name; optionsPrefix scopes its petscOptions
    [[nodiscard]] static std::unique_ptr<LinearSolver> create
    (
        const Name& name,
        Scalar tolerance,
        Count maxIterations,
        const Name& optionsPrefix
    );

    /// Names of every selectable linear solver
    [[nodiscard]] static NameList availableSolvers();

// ***************************** Accessor Methods *****************************

    /// Get relative tolerance
    [[nodiscard]] Scalar tolerance() const noexcept
    {
        return tolerance_;
    }

    /// Get maximum iterations
    [[nodiscard]] Count maxIterations() const noexcept
    {
        return maxIterations_;
    }

    /// Full diagnostics from the last solve call
    [[nodiscard]] const SolvePerformance& lastPerformance() const noexcept
    {
        return lastPerformance_;
    }

// ******************************* Solver Method ******************************

    /// Solve the sparse system, writing the result in place through x
    virtual void solve
    (
        std::span<Scalar> x,
        Mat A,
        Vec b
    ) = 0;

    /// Solver name, used in diagnostic output
    [[nodiscard]] const Name& name() const noexcept
    {
        return name_;
    }

// ***************************** Protected Methods ****************************

protected:

    /// Store diagnostics from the most recent solve call
    void setLastPerformance(const SolvePerformance& performance) noexcept
    {
        lastPerformance_ = performance;
    }

// ****************************** Private Members *****************************

private:

    /// Solver name
    Name name_;

    /// Relative residual tolerance for convergence
    Scalar tolerance_;

    /// Maximum solver iterations before failure
    Count maxIterations_;

    /// Diagnostics from the most recent solve call
    SolvePerformance lastPerformance_
    {
        .solverName    = {},
        .iterations    = -1,
        .finalResidual = std::numeric_limits<Scalar>::quiet_NaN(),
        .converged     = false
    };
};


// ************************** class PetscLinearSolver *************************

class PetscLinearSolver final : public LinearSolver
{
public:

    /// Construct the owned KSP with the given Krylov type
    PetscLinearSolver
    (
        Name name,
        KSPType kspType,
        Scalar tolerance,
        Count maxIterations,
        const Name& optionsPrefix
    );

    /// Destructor releases the KSP and the solution wrapper
    ~PetscLinearSolver() noexcept override;

    /// Solve the sparse system using the owned KSP
    void solve
    (
        std::span<Scalar> x,
        Mat A,
        Vec b
    ) final;

// ****************************** Private Members *****************************

private:

    /// Owned PETSc Krylov solver
    KSP ksp_ = nullptr;

    /// Reusable solution vector, bound over the caller's storage per solve
    Vec solution_ = nullptr;

    /// Previous iterate for the non-finite rollback guard
    ScalarList previousSolution_;
};
