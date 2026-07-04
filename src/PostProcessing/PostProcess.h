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
#include "Scalar.h"
#include "StringTypes.h"

// *************************** Forward Declarations ***************************

class TurbulenceModel;
class Mesh;
class MomentumTransport;
struct CaseConfiguration;

namespace VTK
{
class HDF5CellData;
class HDF5BoundaryData;
}

// *************************** namespace PostProcess **************************

namespace PostProcess
{

/// Extract solution fields and print flow statistics
void reportStatistics(const MomentumTransport& solver);

/// Write the CellData and BoundaryData VTKHDF results
void exportResults
(
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
);

/// CellData .vtkhdf path derived from the output file name
[[nodiscard]] FilePath cellDataPath(const CaseConfiguration& config);

/// BoundaryData _boundary.vtkhdf path derived from the output file name
[[nodiscard]] FilePath boundaryDataPath(const CaseConfiguration& config);

/// Gather the solution fields and append one time step to both writers
void appendTimeStep
(
    VTK::HDF5CellData& volumeWriter,
    VTK::HDF5BoundaryData& boundaryWriter,
    Scalar time,
    const MomentumTransport& solver,
    const TurbulenceModel& turbulence
);

} // namespace PostProcess
