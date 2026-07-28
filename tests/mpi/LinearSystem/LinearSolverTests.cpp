/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LinearSolverTests.cpp
 * @brief Both Krylov solvers reproduce a known solution on a tiny SPD system
 *
 * @details A 1D pure-diffusion problem on an 8x2x2 box decomposed across the
 * ranks. The exact solution is the linear profile phi(x) = x, which
 * finite-volume diffusion reproduces at the cell centres. The system is
 * assembled once and solved with each Krylov solver, so the profile must come
 * out the same for either solver at any rank count.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <span>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "Matrix.h"
#include "TransportEquation.h"
#include "LinearSolvers.h"
#include "BoundaryConditions.h"
#include "FixedValue.h"
#include "ZeroGradient.h"
#include "LeastSquares.h"
#include "MeshFixtures.h"
#include "CellData.h"
#include "FaceData.h"
#include "Field.h"
#include "Comm.h"
#include "StringTypes.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;

// ***************************** Internal Helpers *****************************

namespace
{

/// Register the fixed-value / zero-gradient BCs of the 1D diffusion problem
void registerDiffusionBoundaries(BoundaryConditions& bc, Mesh& mesh)
{
    for (const BoundaryPatch& patch : mesh.patches())
    {
        bc.addPatch(patch);
    }

    bc.linkFaces(mesh.faces());

    bc.setBoundaryType
    (
        BoxPatch::xMin,
        Field::p,
        std::make_unique<FixedValue>(Field::p, S(0.0))
    );
    bc.setBoundaryType
    (
        BoxPatch::xMax,
        Field::p,
        std::make_unique<FixedValue>(Field::p, S(8.0))
    );

    for
    (
        const Name& lateral :
        {
            BoxPatch::yMin,
            BoxPatch::yMax,
            BoxPatch::zMin,
            BoxPatch::zMax
        }
    )
    {
        bc.setBoundaryType
        (
            lateral,
            Field::p,
            std::make_unique<ZeroGradient>(Field::p)
        );
    }

    bc.finalize();
}

} // namespace

// *********************** Krylov Solvers On A Box **************************

TEST_CASE("Both Krylov solvers reproduce the linear profile", "[petsc]")
{
    DecomposedBoxMesh box(8, 2, 2);

    BoundaryConditions bc;
    registerDiffusionBoundaries(bc, box.mesh());

    const LeastSquares gradScheme(box.mesh(), bc);

    ScalarField phi;
    const FaceFluxField gammaFace(S(1.0));
    const ScalarField source;
    const VectorField gradPhi;

    const TransportEquation equation
    {
        .field = Field::p,
        .phi = phi,
        .GammaFace = gammaFace,
        .source = source,
        .gradPhi = gradPhi,
        .gradScheme = gradScheme
    };

    Matrix matrix(box.mesh(), bc);
    matrix.buildMatrix(equation);
    matrix.assemble();

    for (const Name& solverName : {Name("PCG"), Name("BiCGSTAB")})
    {
        // A fresh zero-initialised solution vector per solver
        ScalarField solution;
        std::span<Scalar> x(solution.data(), box.mesh().numOwnedCells());

        const auto solver = LinearSolver::create
        (
            solverName, TestTolerances::solverTolerance, Count{200}, "test"
        );
        solver->solve(x, matrix.matrixA(), matrix.rhsVec());

        // CHECK, not REQUIRE: the next iteration creates and runs another
        // Krylov solver, and unwinding one rank out of the loop would leave
        // the others waiting inside those collectives
        CHECK(solver->lastPerformance().converged);
        CHECK(solver->lastPerformance().solverName == solverName);

        for
        (
            Index cellIdx = 0;
            cellIdx < box.mesh().numOwnedCells();
            ++cellIdx
        )
        {
            CHECK_THAT
            (
                solution[cellIdx],
                WithinAbs
                (
                    box.mesh().cells()[cellIdx].centroid().x(),
                    TestTolerances::absSolve
                )
            );
        }
    }
}