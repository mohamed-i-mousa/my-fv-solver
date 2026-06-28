/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PvdTimeSeries.cpp
 * @brief Implementation of PVD collection file helpers
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PvdTimeSeries.h"

// Standard library headers
#include <fstream>
#include <iostream>

// Project headers
#include "ErrorHandler.h"

// ******************************* namespace VTK ******************************

namespace VTK
{

void writePVDTimeSeriesHeader
(
    const FilePath& pvdFile
)
{
    std::ofstream pvdFileOutput(pvdFile);

    if (!pvdFileOutput.is_open())
    {
        FatalError("Failed to open PVD file: " + pvdFile);
    }

    pvdFileOutput
        << "<?xml version=\"1.0\"?>" << '\n'
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << '\n'
        << "  <Collection>" << '\n';

    pvdFileOutput.close();

    std::cout
        << "PVD time series header written: "
        << pvdFile << '\n';
}


void appendPVDTimeStep
(
    const FilePath& pvdFile,
    const FilePath& dataFile,
    Scalar timeValue,
    Count part
)
{
    // Append a single DataSet line to the open collection, O(1), no read-back
    std::ofstream pvdFileOutput(pvdFile, std::ios::app);
    if (!pvdFileOutput.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for appending: "
          + pvdFile
        );
    }

    pvdFileOutput
        << "    <DataSet timestep=\"" << timeValue
        << "\" part=\"" << part
        << "\" file=\"" << dataFile << "\"/>" << '\n';

    pvdFileOutput.close();

    if (pvdFileOutput.fail())
    {
        FatalError("Failed to write PVD file: " + pvdFile);
    }
}


void closePVDTimeSeries(const FilePath& pvdFile)
{
    std::ofstream pvdFileOutput(pvdFile, std::ios::app);
    if (!pvdFileOutput.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for finalizing: "
          + pvdFile
        );
    }

    pvdFileOutput
        << "  </Collection>" << '\n'
        << "</VTKFile>" << '\n';

    pvdFileOutput.close();

    if (pvdFileOutput.fail())
    {
        FatalError("Failed to finalize PVD file: " + pvdFile);
    }
}

} // namespace VTK
