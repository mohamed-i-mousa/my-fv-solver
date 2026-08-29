/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CFDApplication.h
 * @brief Top-level application orchestrator for the CFD solver
 *
 * @details Coordinates the simulation phases in run(). Parsing, setup,
 * solver assembly, reporting, and export.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "StringTypes.h"

// ************************* namespace CFDApplication *************************

namespace CFDApplication
{

/// Run the full simulation
void run(const FilePath& caseFile);

} // namespace CFDApplication
