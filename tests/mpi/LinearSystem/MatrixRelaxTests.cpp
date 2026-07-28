/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MatrixRelaxTests.cpp
 * @brief Patankar under-relaxation algebra in Matrix::relax
 *
 * @details The relaxation algebra is isolated. It is row-local, so every rank
 * stages and checks its own first owned cell and the result is rank-count
 * independent.
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "Matrix.h"
#include "BoundaryConditions.h"
#include "MeshFixtures.h"
#include "CellData.h"
#include "Comm.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ************************* Patankar Under-Relaxation ************************

TEST_CASE("Matrix::relax applies the Patankar algebra", "[petsc]")
{
    DecomposedBoxMesh box(8, 1, 1);
    const BoundaryConditions bc;

    Matrix matrix(box.mesh(), bc);

    // Stage a known diagonal and right-hand side on cell 0
    matrix.diagonal()[0] = S(10.0);
    matrix.vectorB()[0] = S(0.0);

    ScalarField phiPrev;
    phiPrev[0] = S(5.0);

    matrix.relax(S(0.7), phiPrev);

    // a_P -> a_P / alpha
    REQUIRE_THAT
    (
        matrix.diagonal()[0],
        WithinRel(S(10.0) / S(0.7), TestTolerances::relTight)
    );

    // b -> b + ((1 - alpha) / alpha) * a_P * phiPrev
    REQUIRE_THAT
    (
        matrix.vectorB()[0],
        WithinRel
        (
            (S(1.0) - S(0.7)) / S(0.7) * S(10.0) * S(5.0),
            TestTolerances::relTight
        )
    );
}

// ******************************** Identity **********************************

TEST_CASE("Matrix::relax with alpha = 1 is the identity", "[petsc]")
{
    DecomposedBoxMesh box(8, 1, 1);
    const BoundaryConditions bc;

    Matrix matrix(box.mesh(), bc);

    matrix.diagonal()[0] = S(10.0);
    matrix.vectorB()[0] = S(3.0);

    ScalarField phiPrev;
    phiPrev[0] = S(5.0);

    matrix.relax(S(1.0), phiPrev);

    REQUIRE_THAT
    (
        matrix.diagonal()[0],
        WithinRel(S(10.0), TestTolerances::relTight)
    );
    REQUIRE_THAT
    (
        matrix.vectorB()[0],
        WithinRel(S(3.0), TestTolerances::relTight)
    );
}