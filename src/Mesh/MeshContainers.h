/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshContainers.h
 * @brief Intent-revealing aliases for mesh collections
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>

// *************************** Forward Declarations ***************************

class Vector;
class Face;
class Cell;
class BoundaryPatch;

// ********************************** Aliases *********************************

/// An ordered list of node coordinates
using NodeList = std::vector<Vector>;

/// An ordered list of faces
using FaceList = std::vector<Face>;

/// An ordered list of cells
using CellList = std::vector<Cell>;

/// An ordered list of boundary patches
using PatchList = std::vector<BoundaryPatch>;
