/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file PostProcess.cpp
 * @brief After-solve reporting and result export
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "PostProcess.h"

// Standard library headers
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

// Project headers
#include "CaseConfiguration.h"
#include "DerivedFields.h"
#include "ErrorHandler.h"
#include "Logger.h"
#include "Mesh.h"
#include "SIMPLE.h"
#include "VTK/PvdTimeSeries.h"
#include "VTK/VtkBoundaryWriter.h"
#include "VTK/VtkWriter.h"
#include "TurbulenceModel.h"

// *************************** namespace PostProcess **************************

namespace PostProcess
{

    using ScalarFieldMap = std::map<Name, const ScalarField*>;
    using VectorFieldMap = std::map<Name, std::array<const ScalarField*, 3>>;
    using BoundaryScalarFieldMap = std::map<Name, const FaceData<Scalar>*>;

// ***************************** Internal Helpers *****************************

namespace
{

// Output base path with any trailing ".vtu" extension removed
FilePath outputBase(const FilePath& filename)
{
    FilePath base = filename;
    static constexpr FilePathRef extension = ".vtu";
    if (base.ends_with(extension))
    {
        base.resize(base.size() - extension.size());
    }
    return base;
}

// Text after the last path separator (PVD references are directory-relative)
FilePath fileName(const FilePath& path)
{
    const Index slash = path.find_last_of('/');
    return slash == FilePath::npos ? path : path.substr(slash + 1);
}

// Zero-padded step suffix, e.g. step 42 -> "000042"
Name stepSuffix(Count step)
{
    std::ostringstream suffix;
    suffix << std::setw(6) << std::setfill('0') << step;
    return suffix.str();
}

// Collect the volume and boundary fields and write the .vtu and .vtp files
void writeFields
(
    const FilePath& vtuFilename,
    const FilePath& vtpFilename,
    const Mesh& mesh,
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    bool debug
)
{
    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();

    const ScalarField velocityMag = VTK::velocityMagnitude(Ux, Uy, Uz);

    ScalarFieldMap scalarFieldsToVtk;
    scalarFieldsToVtk["pressure"] = &pressure;
    scalarFieldsToVtk["velocityMagnitude"] = &velocityMag;

    for (const auto& output : turbulence.cellDataOutputs())
    {
        if (output.second != nullptr)
        {
            scalarFieldsToVtk[Name{output.first}] = output.second;
        }
    }

    VectorFieldMap vectorFieldsToVtk;
    vectorFieldsToVtk["velocity"] = {&Ux, &Uy, &Uz};

    VTK::writeVtkUnstructuredGrid
    (
        vtuFilename,
        mesh,
        scalarFieldsToVtk,
        vectorFieldsToVtk,
        debug
    );

    BoundaryScalarFieldMap boundaryScalarFields;

    for (const auto& output : turbulence.boundaryDataOutputs())
    {
        if (output.second != nullptr)
        {
            boundaryScalarFields[Name{output.first}] = output.second;
        }
    }

    const FaceData<Scalar> wallShearStress =
        turbulence.wallShearStress(Ux, Uy, Uz);
    boundaryScalarFields["wallShearStress"] = &wallShearStress;

    VTK::writeBoundaryData
    (
        vtpFilename,
        mesh,
        boundaryScalarFields,
        debug
    );
}

// Derive the "<base>_boundary.vtp" path from a ".vtu" path
FilePath boundaryPath(const FilePath& vtuFilename)
{
    FilePath vtpFilename = vtuFilename;
    const Index dotPos = vtpFilename.rfind(".vtu");

    if (dotPos != FilePath::npos)
    {
        vtpFilename.replace(dotPos, 4, "_boundary.vtp");
    }
    else
    {
        vtpFilename += "_boundary.vtp";
    }

    return vtpFilename;
}

} // namespace

void reportStatistics(const SIMPLE& solver)
{
    std::cout << '\n';
    Logger::sectionHeader("Post-Processing Results");

    const ScalarField& Ux = solver.Ux();
    const ScalarField& Uy = solver.Uy();
    const ScalarField& Uz = solver.Uz();
    const ScalarField& pressure = solver.pressure();

    const ScalarField velocityMag = VTK::velocityMagnitude(Ux, Uy, Uz);

    if (Ux.empty())
    {
        Warning("Solution fields are empty. Skipping statistics.");
        return;
    }

    Scalar maximumVelocity = S(0.0);
    Scalar averageVelocity = S(0.0);
    Scalar maximumPressure = pressure[0];
    Scalar minimumPressure = pressure[0];

    #pragma omp parallel for schedule(static) \
        reduction(max:maximumVelocity, maximumPressure) \
        reduction(min:minimumPressure) \
        reduction(+:averageVelocity)
    for (Index cellIdx = 0; cellIdx < Ux.size(); ++cellIdx)
    {
        const Scalar vmag = velocityMag[cellIdx];
        maximumVelocity = std::max(maximumVelocity, vmag);
        averageVelocity += vmag;

        maximumPressure = std::max(maximumPressure, pressure[cellIdx]);
        minimumPressure = std::min(minimumPressure, pressure[cellIdx]);
    }
    averageVelocity /= S(Ux.size());

    Logger::subsection("Flow statistics");
    Logger::keyValue("Max velocity", maximumVelocity, "m/s");
    Logger::keyValue("Average velocity", averageVelocity, "m/s");
    Logger::keyValue("Pressure min", minimumPressure, "Pa");
    Logger::keyValue("Pressure max", maximumPressure, "Pa");
    Logger::iterationFooter();
}


void exportResults
(
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    const Mesh& mesh,
    const CaseConfiguration& config
)
{
    std::cout << '\n';
    Logger::sectionHeader("Exporting Results");

    FilePath vtuFilename = outputBase(config.vtkOutputFilename) + ".vtu";
    const FilePath vtpFilename = boundaryPath(vtuFilename);

    if (config.debug)
    {
        std::cout
            << '\n' << "Exporting results to VTK UnstructuredGrid..." << '\n';
    }

    writeFields
    (
        vtuFilename,
        vtpFilename,
        mesh,
        solver,
        turbulence,
        config.debug
    );

    Logger::keyValue("Volume field", vtuFilename);
    Logger::keyValue("Boundary field", vtpFilename);
    Logger::iterationFooter();
}


FilePath pvdPathFor(const CaseConfiguration& config)
{
    return outputBase(config.vtkOutputFilename) + ".pvd";
}


void exportTimeStep
(
    const FilePath& pvdFile,
    Scalar time,
    Count step,
    const Mesh& mesh,
    const SIMPLE& solver,
    const TurbulenceModel& turbulence,
    const CaseConfiguration& config
)
{
    const FilePath base = outputBase(config.vtkOutputFilename);
    const FilePath vtuFilename = base + "_" + stepSuffix(step) + ".vtu";
    const FilePath vtpFilename = boundaryPath(vtuFilename);

    writeFields
    (
        vtuFilename,
        vtpFilename,
        mesh,
        solver,
        turbulence,
        config.debug
    );

    // PVD references are relative to the collection file's directory.
    // Volume (.vtu) and wall-boundary (.vtp) share the timestep as parts 0/1
    // so ParaView animates them together as one multiblock.
    VTK::appendPVDTimeStep(pvdFile, fileName(vtuFilename), time, 0);
    VTK::appendPVDTimeStep(pvdFile, fileName(vtpFilename), time, 1);
}

} // namespace PostProcess
