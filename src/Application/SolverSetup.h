/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file SolverSetup.h
 * @brief Runtime service ownership and SIMPLE solver assembly
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <memory>

// *************************** Forward Declarations ***************************

class BoundaryConditions;
class TimeScheme;
class GradientScheme;
class ConvectionScheme;
class TurbulenceModel;
class LinearSolver;
class Mesh;
class MomentumTransport;
struct CaseConfiguration;

// *************************** struct SolverModules ***************************

struct SolverModules
{
    /// Constructor
    SolverModules();

    /// Copy constructor and assignment - Not copyable (unique_ptr members)
    SolverModules(const SolverModules&) = delete;
    SolverModules& operator=(const SolverModules&) = delete;

    /// Move constructor and assignment - Not movable (borrowed references)
    SolverModules(SolverModules&&) = delete;
    SolverModules& operator=(SolverModules&&) = delete;

    /// Destructor
    ~SolverModules() noexcept;

    /// Time-derivative scheme
    std::unique_ptr<TimeScheme> timeScheme;

    /// Gradient scheme owned for SIMPLE's borrowed reference
    std::unique_ptr<GradientScheme> gradScheme;

    /// Default convection scheme fallback
    std::unique_ptr<ConvectionScheme> defaultConvectionScheme;

    /// Momentum equation convection scheme override
    std::unique_ptr<ConvectionScheme> momentumConvectionScheme;

    /// k equation convection scheme override
    std::unique_ptr<ConvectionScheme> kConvectionScheme;

    /// omega equation convection scheme override
    std::unique_ptr<ConvectionScheme> omegaConvectionScheme;

    /// Momentum linear solver
    std::unique_ptr<LinearSolver> momentumSolver;

    /// Pressure-correction linear solver
    std::unique_ptr<LinearSolver> pressureSolver;

    /// k equation linear solver
    std::unique_ptr<LinearSolver> kSolver;

    /// omega equation linear solver
    std::unique_ptr<LinearSolver> omegaSolver;

    /// Turbulence model
    std::unique_ptr<TurbulenceModel> turbulenceModel;

    /// Momentum transport solver
    std::unique_ptr<MomentumTransport> solver;
};

// *************************** namespace SolverSetup **************************

namespace SolverSetup
{

/// Configure runtime and construct the momentum-transport solver
void configure
(
    SolverModules& modules,
    const Mesh& mesh,
    BoundaryConditions& boundaryConditions,
    const CaseConfiguration& config
);

/// Print the solver setup banner
void logSetup
(
    const SolverModules& modules,
    const CaseConfiguration& config
);

} // namespace SolverSetup
