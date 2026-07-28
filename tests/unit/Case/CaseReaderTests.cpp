/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CaseReaderTests.cpp
 * @brief Case-file parsing: typed lookups, sections, defaults
 *
 * @details The reader is exercised against a committed fixture case file
 * (tests/fixtures/cases/parserCase). The fixture directory is pinned at
 * configure time through TURBLYZE_TEST_FIXTURE_DIR, so the test never depends
 * on the working directory ctest happens to launch from.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>
#include <string>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "CaseReader.h"
#include "Vector.h"
#include "Scalar.h"
#include "Integer.h"
#include "StringTypes.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Absolute path of the committed parser fixture case file
[[nodiscard]] FilePath parserCasePath()
{
    return FilePath(TURBLYZE_TEST_FIXTURE_DIR) + "/cases/parserCase";
}

/// Whether a name list carries a given token
[[nodiscard]] bool contains(const NameList& names, const Name& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

} // namespace

// ***************************** Typed Lookups *******************************

TEST_CASE("CaseReader converts each supported value type", "[case]")
{
    const CaseReader reader(parserCasePath());
    const CaseReader& scalars = reader.section("scalars");

    // Scalar, including scientific notation
    REQUIRE_THAT
    (
        scalars.lookup<Scalar>("rho"),
        WithinRel(S(1.225), TestTolerances::relTight)
    );
    REQUIRE_THAT
    (
        scalars.lookup<Scalar>("mu"),
        WithinRel(S(1.7894e-5), TestTolerances::relTight)
    );

    // Count and int (integers assert exactly)
    REQUIRE(scalars.lookup<Count>("iterations") == Count{42});
    REQUIRE(scalars.lookup<int>("offset") == -3);

    // Boolean true- and false-families
    REQUIRE(scalars.lookup<bool>("enabled"));
    REQUIRE_FALSE(scalars.lookup<bool>("verbose"));

    // Token passes through verbatim
    REQUIRE(scalars.lookup<Token>("model") == "kOmegaSST");

    // Vector parses from (x y z)
    const Vector velocity = scalars.lookup<Vector>("velocity");
    REQUIRE_THAT(velocity.x(), WithinRel(S(1.0), TestTolerances::relTight));
    REQUIRE_THAT(velocity.y(), WithinRel(S(2.0), TestTolerances::relTight));
    REQUIRE_THAT(velocity.z(), WithinRel(S(3.0), TestTolerances::relTight));
}

// **************************** Multi-Word Values *****************************

TEST_CASE("CaseReader keeps the separators inside a value", "[case]")
{
    const CaseReader reader(parserCasePath());
    const CaseReader& multiWord = reader.section("multiWord");

    // A path whose directory carries a space stays one usable path
    REQUIRE
    (
        multiWord.lookup<FilePath>("spacedPath")
     == "/tmp/Mobile Documents/sphere.msh"
    );

    // The documented petscOptions form: every token keeps its separator
    REQUIRE
    (
        multiWord.lookup<Token>("options")
     == "-pressure_pc_type icc -momentum_ksp_view"
    );

    // A run of whitespace collapses to exactly one space
    REQUIRE(multiWord.lookup<Token>("collapsed") == "one two");

    // A comment before the ';' is a separator, not part of the value
    REQUIRE(multiWord.lookup<Token>("commented") == "value");
    REQUIRE(multiWord.lookup<Token>("lineCommented") == "value");

    // Mid-word '//' and '/*' are path characters, not comment openers
    REQUIRE
    (
        multiWord.lookup<FilePath>("doubledSlash")
     == "/tmp/My Meshes//v2/sphere.msh"
    );
    REQUIRE
    (
        multiWord.lookup<Token>("starPath") == "My Meshes/*v2*/sphere.msh"
    );
}

TEST_CASE("CaseReader keeps the character after a vector's ')'", "[case]")
{
    const CaseReader reader(parserCasePath());
    const CaseReader& edge = reader.section("vectorEdge");

    const Vector adjacent = edge.lookup<Vector>("adjacent");
    REQUIRE_THAT(adjacent.x(), WithinRel(S(4.0), TestTolerances::relTight));
    REQUIRE_THAT(adjacent.y(), WithinRel(S(5.0), TestTolerances::relTight));
    REQUIRE_THAT(adjacent.z(), WithinRel(S(6.0), TestTolerances::relTight));

    // The key right after the ')' is not truncated
    REQUIRE(edge.lookup<Count>("after") == Count{9});
}

// ****************************** Default Values ******************************

TEST_CASE("CaseReader lookupOrDefault falls back on a missing key", "[case]")
{
    const CaseReader reader(parserCasePath());
    const CaseReader& scalars = reader.section("scalars");

    // Present key ignores the default
    REQUIRE
    (
        scalars.lookupOrDefault<Count>("iterations", Count{7}) == Count{42}
    );

    // Absent key returns the default
    REQUIRE_THAT
    (
        scalars.lookupOrDefault<Scalar>("absent", S(9.5)),
        WithinRel(S(9.5), TestTolerances::relTight)
    );
}

// ***************************** Nested Sections *****************************

TEST_CASE("CaseReader resolves nested sections by chaining", "[case]")
{
    const CaseReader reader(parserCasePath());

    const CaseReader& inner = reader.section("nested").section("inner");

    REQUIRE(inner.lookup<Count>("depth") == Count{2});
    REQUIRE_THAT
    (
        inner.lookup<Scalar>("value"),
        WithinRel(S(7.5), TestTolerances::relTight)
    );
}

// **************************** Section Discovery ****************************

TEST_CASE("CaseReader reports its section names", "[case]")
{
    const CaseReader reader(parserCasePath());

    REQUIRE(reader.hasSection("scalars"));
    REQUIRE(reader.hasSection("boundaryTypes"));
    REQUIRE_FALSE(reader.hasSection("noSuchSection"));

    const NameList names = reader.sectionNames();
    REQUIRE(contains(names, "scalars"));
    REQUIRE(contains(names, "nested"));
    REQUIRE(contains(names, "boundaryTypes"));
}