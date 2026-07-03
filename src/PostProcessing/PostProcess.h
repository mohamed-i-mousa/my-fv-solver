/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PostProcess.h
 * @brief After-solve reporting and result export
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "Integer.h"
#include "Scalar.h"
#include "StringTypes.h"

// *************************** Forward Declarations ***************************

class TurbulenceModel;
class Mesh;
class MomentumTransport;
struct CaseConfiguration;

// *************************** namespace PostProcess **************************

namespace PostProcess
{

/// Extract solution fields and print flow statistics
void reportStatistics(const MomentumTransport& solver);

/// Write VTK volume and wall-boundary results
void exportResults
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
);

/// PVD collection-file path derived from the configured VTK output filename
[[nodiscard]] FilePath pvdPathFor(const CaseConfiguration& config);

/// Write one time step's VTK results and append it to the PVD time series
void exportTimeStep
(
    const FilePath& pvdFile,
    Scalar time,
    Count step,
    const Mesh& mesh,
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const CaseConfiguration& config
);

} // namespace PostProcess
