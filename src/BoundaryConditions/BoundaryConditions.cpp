/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryConditions.cpp
 * @brief Implementation of boundary conditions management system
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "BoundaryConditions.h"

// Standard library headers
#include <iostream>
#include <utility>
#include <set>

// Project headers
#include "ErrorHandler.h"
#include "Reduce.h"

// ***************************** Internal Helpers *****************************

namespace
{

/// Registered boundary condition object for a field, or nullptr
const BoundaryType* findBoundaryType
(
    const std::map<Field, std::unique_ptr<BoundaryType>>& fieldMap,
    Field field
)
{
    const auto fieldIterator = fieldMap.find(field);

    return
        fieldIterator != fieldMap.end()
      ? fieldIterator->second.get()
      : nullptr;
}

} // namespace

// ****************************** Setter Methods ******************************

void BoundaryConditions::addPatch(BoundaryPatch patch)
{
    if (linked_)
    {
        FatalError
        (
            "Cannot add patch after linkFaces() has been called. "
            "Stored face pointers would become invalid."
        );
    }
    patches_.push_back(std::move(patch));
}


void BoundaryConditions::setBoundaryType
(
    const Name& patchName,
    Field field,
    std::unique_ptr<BoundaryType> bc
)
{
    if (finalized_)
    {
        FatalError
        (
            "Cannot register boundary condition for field '"
          + Name(fieldToString(field)) + "' on patch '" + patchName
          + "' after finalize()."
        );
    }

    boundaryTypes_[patchName][field] = std::move(bc);
}

// ***************************** Accessor Methods *****************************

const BoundaryPatch& BoundaryConditions::patch(const Name& patchName) const
{
    for (const auto& patch : patches_)
    {
        if (patch.name() == patchName)
        {
            return patch;
        }
    }

    FatalError("Patch " + patchName + " not found");
}


bool BoundaryConditions::hasBoundaryType
(
    const Name& patchName,
    Field field
) const noexcept
{
    const auto patchIterator = boundaryTypes_.find(patchName);

    return patchIterator != boundaryTypes_.end()
        && patchIterator->second.contains(field);
}


const BoundaryType& BoundaryConditions::boundaryType
(
    const Name& patchName,
    Field field
) const
{
    const auto patchIterator = boundaryTypes_.find(patchName);

    if (patchIterator != boundaryTypes_.end())
    {
        const auto fieldIterator = patchIterator->second.find(field);

        if (fieldIterator != patchIterator->second.end())
        {
            return *fieldIterator->second;
        }
    }

    FatalError
    (
        "Boundary condition not found for patch " + patchName
      + " and field " + Name(fieldToString(field))
    );
}


Index BoundaryConditions::boundaryIdx(Index faceIdx) const
{
    const Index compactIdx = boundaryIdx_[faceIdx];

    if (compactIdx == noBoundaryIdx_)
    {
        FatalError
        (
            "Face " + std::to_string(faceIdx)
          + " is not a physical boundary face."
        );
    }

    return compactIdx;
}


const BoundaryType& BoundaryConditions::boundaryType
(
    Field field,
    Index boundaryIdx
) const
{
    const BoundaryType* type = boundaryTypeAt_[fieldSlot(field)][boundaryIdx];

    if (type == nullptr)
    {
        FatalError
        (
            "No boundary condition registered for field '"
          + Name(fieldToString(field))
          + "' at boundary face " + std::to_string(boundaryIdx) + "."
        );
    }

    return *type;
}


void BoundaryConditions::linkFaces(FaceList& faces)
{
    for (const auto& patch : patches_)
    {
        if (patch.firstFaceIdx() > patch.lastFaceIdx())
        {
            FatalError
            (
                "Boundary patch '" + patch.name()
              + "' has an inverted face range ["
              + std::to_string(patch.firstFaceIdx()) + ", "
              + std::to_string(patch.lastFaceIdx()) + "]."
            );
        }

        if (patch.lastFaceIdx() >= faces.size())
        {
            FatalError
            (
                "Boundary patch '" + patch.name()
              + "' references face index "
              + std::to_string(patch.lastFaceIdx())
              + " outside the valid range [0, "
              + std::to_string(faces.size()) + ")."
            );
        }

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            faces[faceIdx].setPatch(patch);
        }
    }

    linked_ = true;

    // Compact boundary indexing over physical patches; processor patches are
    // inter-rank cuts whose faces carry a neighbor cell, not boundary physics
    boundaryIdx_.assign(faces.size(), noBoundaryIdx_);
    patchStart_.assign(patches_.size(), noBoundaryIdx_);

    Count numBoundaryFaces = 0;

    for (Index patchIdx = 0; patchIdx < patches_.size(); ++patchIdx)
    {
        if (patches_[patchIdx].type() == PatchType::processor)
        {
            continue;
        }

        patchStart_[patchIdx] = numBoundaryFaces;
        numBoundaryFaces += patches_[patchIdx].numFaces();
    }

    geomOwnerCells_.resize(numBoundaryFaces);
    normals_.resize(numBoundaryFaces);
    diffMetric_.resize(numBoundaryFaces);
    normalDistance_.resize(numBoundaryFaces);
    ownerVelocity_.resize(numBoundaryFaces);

    for (Index patchIdx = 0; patchIdx < patches_.size(); ++patchIdx)
    {
        if (patchStart_[patchIdx] == noBoundaryIdx_)
        {
            continue;
        }

        Index compactIdx = patchStart_[patchIdx];

        for
        (
            Index faceIdx = patches_[patchIdx].firstFaceIdx();
            faceIdx <= patches_[patchIdx].lastFaceIdx();
            ++faceIdx, ++compactIdx
        )
        {
            const Face& face = faces[faceIdx];

            if (!face.isBoundary())
            {
                FatalError
                (
                    "Face " + std::to_string(faceIdx)
                  + " on physical patch '"
                  + patches_[patchIdx].name()
                  + "' has a neighbor cell."
                );
            }

            boundaryIdx_[faceIdx] = compactIdx;

            // Over-relaxed orthogonal metric, per unit diffusivity and area
            const Vector Sf = face.normal() * face.projectedArea();
            const Vector ePf = normalized(face.dPf());
            const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePf)) * ePf;

            geomOwnerCells_[compactIdx] = face.ownerCell();
            normals_[compactIdx] = face.normal();
            diffMetric_[compactIdx] =
                magnitude(Ef)
              / (face.projectedArea() * (face.dPfMag() + vSmallValue));
            normalDistance_[compactIdx] = dot(face.dPf(), face.normal());
        }
    }

    for (Index slot = 0; slot < numFields_; ++slot)
    {
        boundaryTypeAt_[slot].assign(numBoundaryFaces, nullptr);
    }

    fluxConstrained_.assign(faces.size(), 0);
    correctsFlux_.assign(faces.size(), 0);
    velocityHullExcluded_.assign(faces.size(), 0);
}


void BoundaryConditions::finalize()
{
    if (!linked_)
    {
        FatalError("BoundaryConditions::finalize called before linkFaces.");
    }

    if (finalized_)
    {
        FatalError("BoundaryConditions::finalize called twice.");
    }

    for (Index patchIdx = 0; patchIdx < patches_.size(); ++patchIdx)
    {
        if (patchStart_[patchIdx] == noBoundaryIdx_)
        {
            continue;
        }

        const BoundaryPatch& patch = patches_[patchIdx];
        const auto patchIterator = boundaryTypes_.find(patch.name());

        if (patchIterator == boundaryTypes_.end())
        {
            continue;
        }

        const auto& fieldMap = patchIterator->second;

        // Record the type pointer for every registered (field, face) pair
        for (const auto& fieldBCPair : fieldMap)
        {
            const Count slot = fieldSlot(fieldBCPair.first);
            std::vector<const BoundaryType*>& fieldTypeAt =
                boundaryTypeAt_[slot];

            for (Index i = 0; i < patch.numFaces(); ++i)
            {
                fieldTypeAt[patchStart_[patchIdx] + i] =
                    fieldBCPair.second.get();
            }
        }

        // Per-face trait flags; the velocity flags demand component agreement
        const BoundaryType* UxType = findBoundaryType(fieldMap, Field::Ux);
        const BoundaryType* UyType = findBoundaryType(fieldMap, Field::Uy);
        const BoundaryType* UzType = findBoundaryType(fieldMap, Field::Uz);

        if ((UxType != nullptr) != (UyType != nullptr)
         || (UxType != nullptr) != (UzType != nullptr)
         || (UxType != nullptr
          && (UxType->constrainsZeroFlux() != UyType->constrainsZeroFlux()
           || UxType->constrainsZeroFlux() != UzType->constrainsZeroFlux()
           || UxType->contributesToLimiterHull()
           != UyType->contributesToLimiterHull()
           || UxType->contributesToLimiterHull()
           != UzType->contributesToLimiterHull())))
        {
            FatalError
            (
                "Velocity component boundary conditions disagree on patch '"
              + patch.name() + "'."
            );
        }

        const BoundaryType* pType = findBoundaryType(fieldMap, Field::p);

        const bool zeroFlux = UxType != nullptr && UxType->constrainsZeroFlux();
        const bool hullExcluded =
            UxType != nullptr && !UxType->contributesToLimiterHull();
        const bool fluxCorrected =
            pType != nullptr && pType->correctsBoundaryFlux();

        for
        (
            Index faceIdx = patch.firstFaceIdx();
            faceIdx <= patch.lastFaceIdx();
            ++faceIdx
        )
        {
            fluxConstrained_[faceIdx] = zeroFlux ? 1 : 0;
            velocityHullExcluded_[faceIdx] = hullExcluded ? 1 : 0;
            correctsFlux_[faceIdx] = fluxCorrected ? 1 : 0;
        }
    }

    finalized_ = true;
}


void BoundaryConditions::snapshotBoundaryVelocity
(
    const ScalarField& Ux,
    const ScalarField& Uy,
    const ScalarField& Uz
)
{
    if (!finalized_)
    {
        FatalError
        (
            "BoundaryConditions::snapshotBoundaryVelocity called before "
            "finalize."
        );
    }

    // The symmetry mirror needs the owner velocity of the just-solved
    // components; single-component consumers read it back per boundary face
    for (Index c = 0; c < ownerVelocity_.size(); ++c)
    {
        const Index owner = geomOwnerCells_[c];
        ownerVelocity_[c] = Vector{Ux[owner], Uy[owner], Uz[owner]};
    }
}


void BoundaryConditions::validatePatchNames() const
{
    // std::set guarantees uniqueness
    std::set<Name> validNames;

    for (const auto& patch : patches_)
    {
        validNames.insert(patch.name());
    }

    // A missing patch is an error only when NO rank has it
    for (const auto& entry : boundaryTypes_)
    {
        if (globalOr(validNames.contains(entry.first)))
        {
            continue;
        }

        Message validList;
        for (const auto& name : validNames)
        {
            if (!validList.empty())
            {
                validList += ", ";
            }
            validList += "'" + name + "'";
        }

        FatalError
        (
            "Boundary condition patch '"
          + entry.first
          + "' does not match any mesh patch on any rank. "
            "Local patch names: " + validList
        );
    }
}


void BoundaryConditions::printSummary() const
{
    std::cout
        << '\n'
        << "--- Boundary Conditions Setup Summary ---" << '\n';

    if (patches_.empty())
    {
        std::cout
            << "  No mesh patches loaded." << '\n';

        return;
    }

    std::cout
        << "Total Mesh Patches Loaded: " << patches_.size()
        << '\n';

    for (const auto& meshPatch : patches_)
    {
        std::cout
            << "  ------------------------------------"
            << '\n';

        std::cout
            << "  Mesh Patch Name         : "
            << meshPatch.name() << '\n';

        std::cout
            << "  Zone ID                 : "
            << meshPatch.zoneIdx() << '\n';

        std::cout
            << "  Number of Faces         : "
            << meshPatch.numFaces() << '\n';

        const auto patchIterator =
            boundaryTypes_.find(meshPatch.name());

        if
        (
            patchIterator != boundaryTypes_.end()
         && !patchIterator->second.empty()
        )
        {
            std::cout
                << "  Configured Physical BCs :" << '\n';

            for (const auto& fieldBCPair : patchIterator->second)
            {
                std::cout
                    << "      Field '"
                    << fieldToString(fieldBCPair.first)
                    << "': Type: ";

                fieldBCPair.second->write(std::cout);

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "  ------------------------------------" << '\n';
}
