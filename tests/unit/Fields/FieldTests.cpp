/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file FieldTests.cpp
 * @brief Unit tests for the Field identifier and its string mapping
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "Field.h"

// **************************** fieldToString Mapping *************************

TEST_CASE("fieldToString covers every field", "[fields]")
{
    REQUIRE(fieldToString(Field::Ux) == "Ux");
    REQUIRE(fieldToString(Field::Uy) == "Uy");
    REQUIRE(fieldToString(Field::Uz) == "Uz");
    REQUIRE(fieldToString(Field::p) == "p");
    REQUIRE(fieldToString(Field::pCorr) == "pCorr");
    REQUIRE(fieldToString(Field::k) == "k");
    REQUIRE(fieldToString(Field::omega) == "omega");
    REQUIRE(fieldToString(Field::nut) == "nut");
}