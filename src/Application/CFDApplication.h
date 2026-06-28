/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CFDApplication.h
 * @brief Top-level application orchestrator for the CFD solver
 *
 * @details CFDApplication owns only the case-file path and coordinates the
 * simulation phases in run(). Phase-specific parsing, setup, solver assembly,
 * reporting, and export live in focused modules.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "StringTypes.h"

// *************************** Forward Declarations ***************************

struct SolverModules;
struct CaseConfiguration;
class Mesh;
class BoundaryConditions;

// *************************** class CFDApplication ***************************

class CFDApplication
{
public:

// ************************* Special Member Functions *************************

    /// Constructor for CFDApplication
    explicit CFDApplication(const FilePath& caseFile);

    /// Copy constructor and assignment - Not copyable
    CFDApplication(const CFDApplication&) = delete;
    CFDApplication& operator=(const CFDApplication&) = delete;

    /// Move constructor and assignment - Not movable
    CFDApplication(CFDApplication&&) = delete;
    CFDApplication& operator=(CFDApplication&&) = delete;

    /// Destructor
    ~CFDApplication() noexcept;

// ******************************** Solver Run ********************************

    /// Run the full simulation
    void run();

// ****************************** Private Methods *****************************

private:

    /// Steady-state path: solve once, report, export
    void runSteady
    (
        SolverModules& modules,
        const Mesh& mesh,
        const BoundaryConditions& bcManager,
        const CaseConfiguration& config
    );

    /// Transient path: time-march, writing output at the configured cadence
    void runTransient
    (
        SolverModules& modules,
        const Mesh& mesh,
        const BoundaryConditions& bcManager,
        const CaseConfiguration& config
    );

// ****************************** Private Members *****************************

    /// Path to case file
    FilePath caseFile_;
};
