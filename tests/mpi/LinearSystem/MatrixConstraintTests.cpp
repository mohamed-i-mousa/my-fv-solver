/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MatrixConstraintTests.cpp
 * @brief Matrix::setValues and explicitJacobiUpdate fixed point
 *
 * @details Both operate on the same 1D diffusion assembly whose solution
 * is phi(x) = x. explicitJacobiUpdate started from that exact field must
 * reproduce it (it is a fixed point of A x = b); setValues must pin a cell to
 * a prescribed value that survives the solve. The box is decomposed across
 * the ranks, so both paths cross a cut: the sweep reads ghost values and the
 * pin propagates to the neighbour rank through the exchanged ghost arrays.
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
#include "HaloExchange.h"
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

/// Build the pure-diffusion transport equation for pressure on a box
[[nodiscard]] TransportEquation diffusionEquation
(
    ScalarField& phi,
    const FaceFluxField& gammaFace,
    const ScalarField& source,
    const VectorField& gradPhi,
    const LeastSquares& gradScheme
)
{
    return TransportEquation
    {
        .field = Field::p,
        .phi = phi,
        .GammaFace = gammaFace,
        .source = source,
        .gradPhi = gradPhi,
        .gradScheme = gradScheme
    };
}

} // namespace

// ************************ Explicit Jacobi Fixed Point **********************

TEST_CASE("explicitJacobiUpdate reproduces the exact solution", "[petsc]")
{
    DecomposedBoxMesh box(8, 2, 2);

    BoundaryConditions bc;
    registerDiffusionBoundaries(bc, box.mesh());

    const LeastSquares gradScheme(box.mesh(), bc);
    ScalarField phi;
    const FaceFluxField gammaFace(S(1.0));
    const ScalarField source;
    const VectorField gradPhi;

    Matrix matrix(box.mesh(), bc);
    matrix.buildMatrix
    (
        diffusionEquation(phi, gammaFace, source, gradPhi, gradScheme)
    );

    // Seed the exact linear profile phi(x) = x, then sweep once. The sweep
    // reads the neighbour values, so the ghosts carry the profile too
    ScalarField phiExact;
    for (Index cellIdx = 0; cellIdx < box.mesh().numOwnedCells(); ++cellIdx)
    {
        phiExact[cellIdx] = box.mesh().cells()[cellIdx].centroid().x();
    }

    exchangeHalos<Scalar>(box.mesh(), {&phiExact});

    ScalarField phiNew;
    matrix.explicitJacobiUpdate(phiExact, phiNew);

    for (Index cellIdx = 0; cellIdx < box.mesh().numOwnedCells(); ++cellIdx)
    {
        REQUIRE_THAT
        (
            phiNew[cellIdx],
            WithinAbs(phiExact[cellIdx],
            TestTolerances::absSolve)
        );
    }
}

// ************************ Fixed-Value Cell Constraint **********************

TEST_CASE("setValues pins a cell through the solve", "[petsc]")
{
    DecomposedBoxMesh box(8, 2, 2);

    BoundaryConditions bc;
    registerDiffusionBoundaries(bc, box.mesh());

    const LeastSquares gradScheme(box.mesh(), bc);
    ScalarField phi;
    const FaceFluxField gammaFace(S(1.0));
    const ScalarField source;
    const VectorField gradPhi;

    Matrix matrix(box.mesh(), bc);
    matrix.buildMatrix
    (
        diffusionEquation(phi, gammaFace, source, gradPhi, gradScheme)
    );

    // Hard-constrain this rank's first cell to a value off the linear profile
    const Index fixedCell = 0;
    const Scalar fixedValue = S(99.0);

    const IndexList cellIndices{fixedCell};
    const ScalarList values{fixedValue};
    const ScalarList fractions;         // empty => full constraint (f = 1)

    // The neighbour rank must see the constraint on its ghost copy, exactly
    // as kOmegaSST exchanges its wall-cell fraction and value
    ScalarField ghostFractions;
    ScalarField ghostValues;
    ghostFractions[fixedCell] = S(1.0);
    ghostValues[fixedCell] = fixedValue;

    exchangeHalos<Scalar>(box.mesh(), {&ghostFractions, &ghostValues});

    matrix.setValues
    (
        cellIndices,
        values,
        ghostFractions,
        ghostValues,
        fractions
    );
    matrix.assemble();

    ScalarField solution;
    std::span<Scalar> x(solution.data(), box.mesh().numOwnedCells());

    const auto solver = LinearSolver::create
    (
        "PCG",
        "Jacobi",
        TestTolerances::solverTolerance,
        Count{200},
        "constraint"
    );
    solver->solve(x, matrix.matrixA(), matrix.rhsVec());

    REQUIRE(solver->lastPerformance().converged);
    REQUIRE_THAT
    (
        solution[fixedCell],
        WithinAbs(fixedValue, TestTolerances::absSolve)
    );
}