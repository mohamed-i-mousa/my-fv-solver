/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HaloExchangeTests.cpp
 * @brief Ghost cells receive their owning rank's values across a cut
 *
 * @details Each rank owns one cell of a 1D chain whose global index equals
 * the rank. The owned value convention f(id) = 100 + id lets the exchange be
 * checked analytically. After exchangeHalos, ghost cell p must hold
 * 100 + neighborRank(p).
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "HaloExchange.h"
#include "MeshFixtures.h"
#include "CellData.h"
#include "Vector.h"
#include "Comm.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Owned-cell value convention: f(globalIdx) = base + globalIdx
[[nodiscard]] Scalar chainValue(Scalar base, Index globalId) noexcept
{
    return base + S(globalId);
}

} // namespace

// **************************** Scalar Ghost Fill ****************************

TEST_CASE("exchangeHalos fills scalar ghosts from the neighbour", "[mpi][parallel]")
{
    if (!Comm::parallelRun())
    {
        SKIP("halo exchange is not available for a single rank");
    }

    const DecomposedChainMesh chain;
    const Mesh& mesh = chain.mesh();

    ScalarField phi;
    phi[0] = chainValue(S(100.0), Comm::myProcessorNum());

    exchangeHalos<Scalar>(mesh, {&phi});

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        REQUIRE_THAT
        (
            phi[patch.ghostFirstCell()],
            WithinRel
            (
                chainValue(S(100.0),
                patch.neighborRank()),
                TestTolerances::relTight
            )
        );
    }
}

// **************************** Vector Ghost Fill ****************************

TEST_CASE
(
    "exchangeHalos fills vector ghosts from the neighbour",
    "[mpi][parallel]"
)
{
    if (!Comm::parallelRun())
    {
        SKIP("halo exchange is a no-op on a single rank");
    }

    const DecomposedChainMesh chain;
    const Mesh& mesh = chain.mesh();

    VectorField velocity;
    velocity[0] =
        Vector(chainValue(S(100.0), Comm::myProcessorNum()), S(0.0), S(0.0));

    exchangeHalos<Vector>(mesh, {&velocity});

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        const Vector& ghost = velocity[patch.ghostFirstCell()];
        REQUIRE_THAT
        (
            ghost.x(),
            WithinRel
            (
                chainValue(S(100.0),
                patch.neighborRank()),
                TestTolerances::relTight
            )
        );
        REQUIRE(ghost.y() == S(0.0));
        REQUIRE(ghost.z() == S(0.0));
    }
}

// ************************** Batched Two-Field Fill *************************

TEST_CASE("exchangeHalos fills two batched fields at once", "[mpi][parallel]")
{
    if (!Comm::parallelRun())
    {
        SKIP("halo exchange is a no-op on a single rank");
    }

    const DecomposedChainMesh chain;
    const Mesh& mesh = chain.mesh();

    const Index rank = Comm::myProcessorNum();
    ScalarField a;
    ScalarField b;
    a[0] = chainValue(S(100.0), rank);
    b[0] = chainValue(S(200.0), rank);

    // One message per neighbour carries both fields, laid out [a ghosts|b ghosts]
    exchangeHalos<Scalar>(mesh, {&a, &b});

    for (const ProcessorPatch& patch : mesh.processorPatches())
    {
        REQUIRE_THAT
        (
            a[patch.ghostFirstCell()],
            WithinRel
            (
                chainValue(S(100.0),
                patch.neighborRank()),
                TestTolerances::relTight
            )
        );
        REQUIRE_THAT
        (
            b[patch.ghostFirstCell()],
            WithinRel
            (
                chainValue(S(200.0),
                patch.neighborRank()),
                TestTolerances::relTight
            )
        );
    }
}