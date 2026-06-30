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

// ***************************** Internal Helpers *****************************

namespace
{
    const Message pvdClosingBlock =
        "  </Collection>\n"
        "</VTKFile>\n";
}

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
        << "  <Collection>" << '\n'
        << pvdClosingBlock;

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
    // Open for in-place update so the closing tags can be overwritten
    std::fstream pvdFileOutput(pvdFile, std::ios::in | std::ios::out);
    if (!pvdFileOutput.is_open())
    {
        FatalError
        (
            "Failed to open PVD file for appending: "
          + pvdFile
        );
    }

    // Rewind over the closing block, insert the DataSet, re-emit the close
    pvdFileOutput.seekp(0, std::ios::end);
    const std::streamoff fileEnd = pvdFileOutput.tellp();
    pvdFileOutput.seekp
    (
        fileEnd - static_cast<std::streamoff>(pvdClosingBlock.size())
    );

    pvdFileOutput
        << "    <DataSet timestep=\"" << timeValue
        << "\" part=\"" << part
        << "\" file=\"" << dataFile << "\"/>" << '\n'
        << pvdClosingBlock;

    pvdFileOutput.close();

    if (pvdFileOutput.fail())
    {
        FatalError("Failed to write PVD file: " + pvdFile);
    }
}

} // namespace VTK
