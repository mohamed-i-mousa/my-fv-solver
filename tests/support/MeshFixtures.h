/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshFixtures.h
 * @brief Programmatic hex-box meshes for unit tests

 * @details The smallest shipped mesh is far too large to be a unit fixture,
 * so tests build tiny structured hex boxes in memory.
 *
 * RAII owner of one populated Mesh: releases the global cell/face counts on
 * destruction so the next TEST_CASE can build its own mesh.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "Mesh.h"
#include "Scalar.h"
#include "Integer.h"

// ********************************** BoxPatch ********************************

// The six axis-aligned boundary patches of a hex box, in emission order
namespace BoxPatch
{

inline const Name xMin = "xMin";
inline const Name xMax = "xMax";
inline const Name yMin = "yMin";
inline const Name yMax = "yMax";
inline const Name zMin = "zMin";
inline const Name zMax = "zMax";

} // namespace BoxPatch

// **************************** Fixture Factories *****************************

/// Build a structured hex box of nx*ny*nz unit-spacing cells at the origin.
[[nodiscard]] Mesh makeHexBoxMesh
(
    Count nx,
    Count ny,
    Count nz,
    Scalar spacing = S(1.0)
);

/// Build this rank's submesh of a 1D cell chain: one owned cell whose global
/// index equals the rank, plus one ghost stub per neighbouring rank. Only the
/// data exchangeHalos reads (processor patches, owned/ghost cell counts) is
/// populated; nodes, faces, and cell geometry are intentionally empty. Every
/// rank must build it together (the ghost count follows the rank position).
[[nodiscard]] Mesh makeDecomposedChainMesh();

/// Build this rank's submesh of the same hex box, partitioned by METIS through
/// the production decompose-distribute path. Serial runs get the whole box, so
/// np = 1 is bit-identical to makeHexBoxMesh. Collective: every rank must
/// build it together.
[[nodiscard]] Mesh makeDecomposedHexBoxMesh
(
    Count nx,
    Count ny,
    Count nz,
    Scalar spacing = S(1.0)
);

// ****************************** class TestMesh ******************************

class TestMesh
{
public:

// ************************* Special Member Functions *************************

    /// Build a hex box
    TestMesh(Count nx, Count ny, Count nz, Scalar spacing = S(1.0))
    :
        mesh_(makeHexBoxMesh(nx, ny, nz, spacing))
    {}

    /// Not copyable movable
    TestMesh(const TestMesh&) = delete;
    TestMesh& operator=(const TestMesh&) = delete;

    /// Not movable
    TestMesh(TestMesh&&) = delete;
    TestMesh& operator=(TestMesh&&) = delete;

    /// Release the global cell/face counts
    ~TestMesh() noexcept
    {
        Mesh::resetCounts();
    }

// ***************************** Accessor Methods *****************************

    /// The owned mesh
    [[nodiscard]] const Mesh& mesh() const noexcept
    {
        return mesh_;
    }

    /// The owned mesh (mutable for face-patch linking)
    [[nodiscard]] Mesh& mesh() noexcept
    {
        return mesh_;
    }

// ****************************** Private Members *****************************

private:

    /// The single populated mesh
    Mesh mesh_;
};

// ************************* class DecomposedChainMesh ************************

class DecomposedChainMesh
{
public:

// ************************* Special Member Functions *************************

    /// Build this rank's chain submesh
    DecomposedChainMesh()
    :
        mesh_(makeDecomposedChainMesh())
    {}

    /// Not copyable
    DecomposedChainMesh(const DecomposedChainMesh&) = delete;
    DecomposedChainMesh& operator=(const DecomposedChainMesh&) = delete;

    /// Not movable
    DecomposedChainMesh(DecomposedChainMesh&&) = delete;
    DecomposedChainMesh& operator=(DecomposedChainMesh&&) = delete;

    /// Release the global cell/face counts
    ~DecomposedChainMesh() noexcept
    {
        Mesh::resetCounts();
    }

// ***************************** Accessor Methods *****************************

    /// The owned submesh
    [[nodiscard]] const Mesh& mesh() const noexcept
    {
        return mesh_;
    }

// ****************************** Private Members *****************************

private:

    /// The single populated submesh
    Mesh mesh_;
};

// ************************** class DecomposedBoxMesh *************************

class DecomposedBoxMesh
{
public:

// ************************* Special Member Functions *************************

    /// Build this rank's share of a hex box
    DecomposedBoxMesh(Count nx, Count ny, Count nz, Scalar spacing = S(1.0))
    :
        mesh_(makeDecomposedHexBoxMesh(nx, ny, nz, spacing))
    {}

    /// Not copyable
    DecomposedBoxMesh(const DecomposedBoxMesh&) = delete;
    DecomposedBoxMesh& operator=(const DecomposedBoxMesh&) = delete;

    /// Not movable
    DecomposedBoxMesh(DecomposedBoxMesh&&) = delete;
    DecomposedBoxMesh& operator=(DecomposedBoxMesh&&) = delete;

    /// Release the global cell/face counts
    ~DecomposedBoxMesh() noexcept
    {
        Mesh::resetCounts();
    }

// ***************************** Accessor Methods *****************************

    /// The owned submesh
    [[nodiscard]] const Mesh& mesh() const noexcept
    {
        return mesh_;
    }

    /// The owned submesh (mutable for face-patch linking)
    [[nodiscard]] Mesh& mesh() noexcept
    {
        return mesh_;
    }

// ****************************** Private Members *****************************

private:

    /// The single populated submesh
    Mesh mesh_;
};