/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file StringTypes.h
 * @brief Intent-revealing string aliases
 *
 * @details Defines aliases over std::string that signal what the text
 * represents, a name, a parser token, a file path, or display text.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <string>
#include <vector>

// ********************************** Aliases *********************************

/// A name identifier
using Name = std::string;

/// A parser token read from a case or mesh file
using Token = std::string;

/// A path or filename
using FilePath = std::string;

/// Human-readable diagnostic or display text
using Message = std::string;

// *********************************** Lists **********************************

/// An ordered collection of owned name identifiers
using NameList = std::vector<Name>;
