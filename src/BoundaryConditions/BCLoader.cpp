/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BCLoader.cpp
 * @brief Case-file boundary condition registration
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BCLoader.h"

// Standard library headers
#include <iostream>
#include <map>
#include <memory>

// Project headers
#include "CaseReader.h"
#include "ErrorHandler.h"
#include "FixedValue.h"
#include "Logger.h"
#include "Reduce.h"
#include "RuntimeSelection.h"
#include "Symmetry.h"
#include "TurbulenceModel.h"
#include "kOmegaSST.h"

// **************************** namespace BCLoader ****************************

namespace BCLoader
{

// ***************************** Internal Helpers *****************************

namespace
{

[[noreturn]] void unknownTypeToken
(
    const Name& bcType,
    const Name& fieldName,
    const Name& patchName,
    const Message& validList
)
{
    FatalError
    (
        "Unknown boundary condition type '" + bcType
      + "' for field '" + fieldName
      + "' on patch '" + patchName
      + "'. Valid types: " + validList
    );
}


/// Validate the type token against the field's selectable names, then
/// construct and register the boundary condition
void registerBC
(
    BoundaryConditions& bcManager,
    const Name& patchName,
    Field field,
    const Name& bcType,
    const CaseReader& patchBC,
    const Name& sectionName
)
{
    if (!RuntimeSelection::isKnown(bcType, BoundaryType::availableTypes(field)))
    {
        unknownTypeToken
        (
            bcType,
            sectionName,
            patchName,
            RuntimeSelection::joinNames(BoundaryType::availableTypes(field))
        );
    }

    bcManager.setBoundaryType(patchName, field, BoundaryType::create(bcType, field, patchBC));
}


void validateWallFunctionSetup
(
    const Mesh& mesh,
    const BoundaryConditions& bcManager,
    const CaseConfiguration& config
)
{
    if (TurbulenceModel::isLaminar(config.turbulenceModel))
    {
        return;
    }

    for (const auto& patch : mesh.patches())
    {
        if (patch.type() != PatchType::wall)
        {
            continue;
        }

        const Name& patchName = patch.name();

        const bool kIsWF =
            bcManager.boundaryType(patchName, Field::k).isWallModelled();
        const bool omegaIsWF =
            bcManager.boundaryType(patchName, Field::omega).isWallModelled();
        const bool nutIsWF =
            bcManager.boundaryType(patchName, Field::nut).isWallModelled();

        const int wfCount = int(kIsWF) + int(omegaIsWF) + int(nutIsWF);

        if (wfCount == 0 || wfCount == 3)
        {
            continue;
        }

        FatalError
        (
            "Wall patch '" + patchName
          + "': wall functions must be configured as a complete triplet "
            "(k + omega + nut) or omitted entirely. Found: k="
          + (kIsWF     ? "WF" : "non-WF")
          + ", omega=" + (omegaIsWF ? "WF" : "non-WF")
          + ", nut="   + (nutIsWF   ? "WF" : "non-WF") + "."
        );
    }
}


bool isSymmetryPatch(const Mesh& mesh, const Name& patchName)
{
    bool symmetryHere = false;

    for (const auto& patch : mesh.patches())
    {
        if (patch.name() == patchName)
        {
            symmetryHere = patch.type() == PatchType::symmetry;
            break;
        }
    }

    // A patch may live on some ranks only
    return globalOr(symmetryHere);
}


void applySymmetry(const Mesh& mesh, BoundaryConditions& bcManager)
{
    using enum Field;
    static constexpr Field solvedFields[] =
        { Ux, Uy, Uz, p, pCorr, k, omega, nut };

    // A symmetry plane is mesh-derived: the case file carries no entry for it
    for (const auto& patch : mesh.patches())
    {
        if (patch.type() != PatchType::symmetry)
        {
            continue;
        }

        for (const Field field : solvedFields)
        {
            bcManager.setBoundaryType
            (
                patch.name(),
                field,
                std::make_unique<Symmetry>(field)
            );
        }
    }
}

} // namespace (unnamed)


// *********************************** Load ***********************************

void load
(
    const CaseReader& reader,
    const CaseConfiguration& config,
    Mesh& mesh,
    BoundaryConditions& bcManager
)
{
    std::cout << '\n';
    Logger::sectionHeader("Setting Boundary Conditions");

    for (const auto& patch : mesh.patches())
    {
        bcManager.addPatch(patch);
    }

    bcManager.linkFaces(mesh.faces());

    for (const auto& face : mesh.faces())
    {
        if (face.isBoundary() && !face.patch().has_value())
        {
            FatalError
            (
                "Boundary face "
              + std::to_string(face.idx())
              + " has no patch after linking."
            );
        }
    }

    const auto& BCs = reader.section("boundaryConditions");

    if (BCs.hasSection("U"))
    {
        const auto& velocityBCs = BCs.section("U");

        for (const auto& patchName : velocityBCs.sectionNames())
        {
            if (isSymmetryPatch(mesh, patchName))
            {
                continue;
            }

            const auto& patchBC = velocityBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            // The case-file vector entry fans out into the scalar components
            registerBC(bcManager, patchName, Field::Ux, bcType, patchBC, "U");
            registerBC(bcManager, patchName, Field::Uy, bcType, patchBC, "U");
            registerBC(bcManager, patchName, Field::Uz, bcType, patchBC, "U");
        }
    }

    bool hasFixedPressure = false;

    if (BCs.hasSection("p"))
    {
        const auto& pressureBCs = BCs.section("p");

        for (const auto& patchName : pressureBCs.sectionNames())
        {
            if (isSymmetryPatch(mesh, patchName))
            {
                continue;
            }

            const auto& patchBC = pressureBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            registerBC(bcManager, patchName, Field::p, bcType, patchBC, "p");

            const BoundaryType& pType = bcManager.boundaryType(patchName, Field::p);

            hasFixedPressure = hasFixedPressure || pType.fixesValue();

            // Derive the p' boundary condition from p: fixed p becomes
            // p' = 0, zero-gradient p stays zero-gradient
            bcManager.setBoundaryType
            (
                patchName,
                Field::pCorr,
                pType.pressureCorrectionCompanion()
            );
        }
    }

    if (!hasFixedPressure)
    {
        Warning
        (
            "No fixedValue pressure boundary condition found. "
            "The pressure field has no reference value, which "
            "may cause a singular pressure matrix."
        );
    }

    // Inlet k per patch for omega's 'calculated' entries: a fixed value
    // carries itself, anything else falls back to the configured inlet k
    std::map<Name, Scalar> resolvedK;

    if (BCs.hasSection("k"))
    {
        const auto& kBCs = BCs.section("k");

        for (const auto& patchName : kBCs.sectionNames())
        {
            if (isSymmetryPatch(mesh, patchName))
            {
                continue;
            }

            const auto& patchBC = kBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            // 'calculated' resolves loader-side: it needs the case config
            if (bcType == "fixedValue")
            {
                const Token valStr = patchBC.lookup<Token>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    value =
                        kOmegaSST::inletK
                        (
                            config.initialVelocity,
                            config.turbulenceIntensity
                        );

                    std::cout
                        << "Inlet turbulence kinetic energy : " << value
                        << '\n';
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager.setBoundaryType
                (
                    patchName,
                    Field::k,
                    std::make_unique<FixedValue>(Field::k, value)
                );
                resolvedK[patchName] = value;
                continue;
            }

            registerBC(bcManager, patchName, Field::k, bcType, patchBC, "k");
            resolvedK[patchName] =
                kOmegaSST::inletK
                (
                    config.initialVelocity,
                    config.turbulenceIntensity
                );
        }
    }

    if (BCs.hasSection("omega"))
    {
        const auto& omegaBCs = BCs.section("omega");

        for (const auto& patchName : omegaBCs.sectionNames())
        {
            if (isSymmetryPatch(mesh, patchName))
            {
                continue;
            }

            const auto& patchBC = omegaBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            // 'calculated' resolves loader-side from the patch's inlet k
            if (bcType == "fixedValue")
            {
                const Token valStr = patchBC.lookup<Token>("value");

                Scalar value = S(0.0);

                if (valStr == "calculated")
                {
                    const auto kIterator = resolvedK.find(patchName);

                    if (kIterator == resolvedK.end())
                    {
                        FatalError
                        (
                            "Boundary condition not found for patch "
                          + patchName + " and field k"
                        );
                    }

                    value =
                        kOmegaSST::inletOmega
                        (
                            kIterator->second,
                            config.hydraulicDiameter
                        );

                    std::cout
                        << "Inlet specific dissipation : " << value << '\n';
                }
                else
                {
                    value = patchBC.lookup<Scalar>("value");
                }

                bcManager.setBoundaryType
                (
                    patchName,
                    Field::omega,
                    std::make_unique<FixedValue>(Field::omega, value)
                );
                continue;
            }

            registerBC
            (
                bcManager,
                patchName,
                Field::omega,
                bcType,
                patchBC,
                "omega"
            );
        }
    }

    if (BCs.hasSection("nut"))
    {
        const auto& nutBCs = BCs.section("nut");

        for (const auto& patchName : nutBCs.sectionNames())
        {
            if (isSymmetryPatch(mesh, patchName))
            {
                continue;
            }

            const auto& patchBC = nutBCs.section(patchName);
            const Name bcType = patchBC.lookup<Name>("type");

            registerBC(bcManager, patchName, Field::nut, bcType, patchBC, "nut");
        }
    }

    bcManager.validatePatchNames();

    validateWallFunctionSetup(mesh, bcManager, config);

    applySymmetry(mesh, bcManager);

    bcManager.finalize();

    if (config.debug)
    {
        bcManager.printSummary();
    }

    std::cout
        << "Boundary conditions set for "
        << mesh.patches().size() << " patches." << '\n';
}

} // namespace BCLoader
