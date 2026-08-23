/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file RuntimeSelectionTests.cpp
 * @brief Factory positive paths and the advertised name lists
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "GradientScheme.h"
#include "ConvectionScheme.h"
#include "TimeScheme.h"
#include "LinearSolvers.h"
#include "TurbulenceModel.h"
#include "BoundaryConditions.h"
#include "MeshFixtures.h"
#include "StringTypes.h"

// ***************************** Internal Helpers *****************************

namespace
{

/// Whether a name list carries a given token
[[nodiscard]] bool contains(const NameList& names, const Name& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

} // namespace

// **************************** Gradient Schemes ******************************

TEST_CASE("Gradient-scheme selection", "[selection]")
{
    const NameList names = GradientScheme::availableSchemes();
    REQUIRE(names.size() == Count{1});
    REQUIRE(names.front() == "leastSquares");

    const TestMesh box(2, 1, 1);
    const BoundaryConditions bc;
    REQUIRE(GradientScheme::create("leastSquares", box.mesh(), bc) != nullptr);
}

// **************************** Convection Schemes ****************************

TEST_CASE("Convection-scheme selection", "[selection]")
{
    const NameList names = ConvectionScheme::availableSchemes();
    REQUIRE(contains(names, "Upwind"));
    REQUIRE(contains(names, "CentralDifference"));
    REQUIRE(contains(names, "SecondOrderUpwind"));

    for (const Name& name : names)
    {
        REQUIRE(ConvectionScheme::create(name) != nullptr);
    }
}

// ****************************** Time Schemes *******************************

TEST_CASE("Time-scheme selection", "[selection]")
{
    const NameList names = TimeScheme::availableSchemes();
    REQUIRE(contains(names, "steadyState"));
    REQUIRE(contains(names, "implicitEuler"));
    REQUIRE(contains(names, "CrankNicolson"));

    for (const Name& name : names)
    {
        REQUIRE(TimeScheme::create(name) != nullptr);
    }
}

// ***************************** Linear Solvers ******************************

TEST_CASE("Linear-solver selection lists every supported method", "[selection]")
{
    const NameList solvers = LinearSolver::availableSolvers();
    REQUIRE(contains(solvers, "BiCGSTAB"));
    REQUIRE(contains(solvers, "PCG"));
    REQUIRE(contains(solvers, "GMRES"));
    REQUIRE(contains(solvers, "FGMRES"));
    REQUIRE(contains(solvers, "TFQMR"));
    REQUIRE(contains(solvers, "CGS"));
    REQUIRE(contains(solvers, "MINRES"));
    REQUIRE(contains(solvers, "Richardson"));
    REQUIRE(contains(solvers, "Chebyshev"));
    REQUIRE(contains(solvers, "PreOnly"));

    const NameList pcs = LinearSolver::availablePreconditioners();
    REQUIRE(contains(pcs, "Jacobi"));
    REQUIRE(contains(pcs, "None"));
    REQUIRE(contains(pcs, "ILU"));
    REQUIRE(contains(pcs, "ICC"));
    REQUIRE(contains(pcs, "SOR"));
    REQUIRE(contains(pcs, "GaussSeidel"));
    REQUIRE(contains(pcs, "AMG"));
    REQUIRE(contains(pcs, "GAMG"));
    REQUIRE(contains(pcs, "BlockJacobi"));
    REQUIRE(contains(pcs, "ASM"));
    REQUIRE(contains(pcs, "LU"));
    REQUIRE(contains(pcs, "Cholesky"));
}

// **************************** Turbulence Models ****************************

TEST_CASE("Turbulence-model selection", "[selection]")
{
    const NameList names = TurbulenceModel::availableModels();
    REQUIRE(contains(names, "Laminar"));
    REQUIRE(contains(names, "kOmegaSST"));

    REQUIRE(TurbulenceModel::isLaminar("Laminar"));
    REQUIRE_FALSE(TurbulenceModel::isLaminar("kOmegaSST"));
}