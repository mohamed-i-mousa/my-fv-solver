/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Matrix.cpp
 * @brief Matrix assembly and linear system construction for equations
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Matrix.h"

// Standard library headers
#include <algorithm>

// Project headers
#include "ErrorHandler.h"
#include "GlobalIndex.h"
#include "PETScRuntime.h"
#include "TimeScheme.h"

// ************************* Special Member Functions *************************

Matrix::Matrix
(
    const Mesh& mesh,
    const BoundaryConditions& boundaryConds
)
:
    mesh_{mesh},
    bcManager_{boundaryConds}
{
    const Count numCells = mesh_.numCells();
    const Count numFaces = mesh_.numFaces();

    // Slot-ordered face lists; faceSlot_ maps a face to its first
    // off-diagonal slot (internal) or its boundary scratch ordinal
    faceSlot_.assign(numFaces, 0);

    for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
    {
        if (mesh_.faces()[faceIdx].isBoundary())
        {
            faceSlot_[faceIdx] = boundaryFaces_.size();
            boundaryFaces_.push_back(faceIdx);
        }
        else
        {
            faceSlot_[faceIdx] = 2 * internalFaces_.size();
            internalFaces_.push_back(faceIdx);
        }
    }

    const Count numInternalFaces = internalFaces_.size();
    const Count numBoundaryFaces = boundaryFaces_.size();

    // COO layout: [off-diagonals | (future processor slots) | diagonals]
    diagOffset_ = 2 * numInternalFaces;
    const Count numCoo = diagOffset_ + numCells;

    // The fixed sparsity pattern: slot 2k couples (owner, neighbor) of
    // internal face k, slot 2k+1 the reverse, and one diagonal per cell.
    // Indices are global through the rank-major cell numbering (identity
    // on one rank), so the pattern is already distribution-ready
    const GlobalIndex globalCells(numCells);

    std::vector<PetscInt> cooRows(numCoo);
    std::vector<PetscInt> cooCols(numCoo);

    for (Index k = 0; k < numInternalFaces; ++k)
    {
        const Face& face = mesh_.faces()[internalFaces_[k]];
        const auto owner =
            static_cast<PetscInt>(globalCells.toGlobal(face.ownerCell()));
        const auto neighbor =
            static_cast<PetscInt>
            (
                globalCells.toGlobal(face.neighborCell().value())
            );

        cooRows[2 * k] = owner;
        cooCols[2 * k] = neighbor;
        cooRows[2 * k + 1] = neighbor;
        cooCols[2 * k + 1] = owner;
    }

    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const auto diagIdx =
            static_cast<PetscInt>(globalCells.toGlobal(cellIdx));

        cooRows[diagOffset_ + cellIdx] = diagIdx;
        cooCols[diagOffset_ + cellIdx] = diagIdx;
    }

    const auto n = static_cast<PetscInt>(numCells);

    PETSC_CHECK(MatCreate(PETScRuntime::comm(), &matrixA_));
    PETSC_CHECK(MatSetSizes(matrixA_, n, n, PETSC_DECIDE, PETSC_DECIDE));
    PETSC_CHECK(MatSetType(matrixA_, MATAIJ));
    PETSC_CHECK
    (
        MatSetPreallocationCOO
        (
            matrixA_,
            static_cast<PetscCount>(numCoo),
            cooRows.data(),
            cooCols.data()
        )
    );

    // Staging and scratch buffers, sized once
    cooValues_.assign(numCoo, S(0.0));
    vectorB_.assign(numCells, S(0.0));
    faceDiag_.assign(2 * numInternalFaces, S(0.0));
    faceRhs_.assign(2 * numInternalFaces, S(0.0));
    boundaryDiag_.assign(numBoundaryFaces, S(0.0));
    boundaryRhs_.assign(numBoundaryFaces, S(0.0));

    // RHS handle wraps vectorB_ storage: staging writes are the vector,
    // with no PETSc-owned backing array allocated behind the wrapper
    // (single-rank form; distributed sizing replaces it in Phase 5)
    PETSC_CHECK
    (
        VecCreateSeqWithArray
        (
            PETScRuntime::comm(),
            1,
            n,
            vectorB_.data(),
            &rhsVec_
        )
    );
}


Matrix::~Matrix() noexcept
{
    if (VecDestroy(&rhsVec_) != PETSC_SUCCESS)
    {
        Warning("VecDestroy failed");
    }
    if (MatDestroy(&matrixA_) != PETSC_SUCCESS)
    {
        Warning("MatDestroy failed");
    }
}

// ****************************** Public Methods ******************************

void Matrix::buildMatrix(const TransportEquation& equation)
{
    lastRelaxationFactor_ = S(0.0);

    const Count numCells = mesh_.numCells();
    const Count numInternalFaces = internalFaces_.size();
    const Count numBoundaryFaces = boundaryFaces_.size();

    // Pass 1: internal faces write their exclusive off-diagonal slots and
    // stage per-face diagonal/RHS contributions for the cell gather
    #pragma omp parallel for schedule(static)
    for (Index k = 0; k < numInternalFaces; ++k)
    {
        assembleInternalFace(k, mesh_.faces()[internalFaces_[k]], equation);
    }

    // Pass 2: boundary faces stage their owner-cell contributions
    #pragma omp parallel for schedule(static)
    for (Index m = 0; m < numBoundaryFaces; ++m)
    {
        assembleBoundaryFace(m, mesh_.faces()[boundaryFaces_[m]], equation);
    }

    // Pass 3: per-cell gather sums the face scratch into the single
    // diagonal slot and RHS entry (race-free: one writer per cell)
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Cell& cell = mesh_.cells()[cellIdx];

        Scalar diag = S(0.0);
        Scalar rhs = equation.source[cellIdx];

        for (const Index faceIdx : cell.faceIndices())
        {
            const Index slot = faceSlot_[faceIdx];

            if (mesh_.faces()[faceIdx].isBoundary())
            {
                diag += boundaryDiag_[slot];
                rhs += boundaryRhs_[slot];
            }
            else
            {
                const Index side =
                    mesh_.faces()[faceIdx].ownerCell() == cellIdx ? 0 : 1;

                diag += faceDiag_[slot + side];
                rhs += faceRhs_[slot + side];
            }
        }

        // Transient term: add d(phi)/dt to the diagonal and RHS
        if (equation.transient)
        {
            const TransientTerm& transient = *equation.transient;

            const Scalar ddtPrevStep =
                transient.ddtPrevStep != nullptr
              ? (*transient.ddtPrevStep)[cellIdx]
              : S(0.0);

            const TimeContribution contribution =
                transient.scheme.coefficients
                (
                    mesh_.cells()[cellIdx].volume(),
                    transient.deltaT,
                    transient.phiPrevStep[cellIdx],
                    ddtPrevStep
                );

            diag += contribution.diag;
            rhs += contribution.source;
        }

        cooValues_[diagOffset_ + cellIdx] = diag;
        vectorB_[cellIdx] = rhs;
    }
}


void Matrix::relax(Scalar alpha, const ScalarField& phiPrevIter)
{
    if (alpha <= S(0.0))
    {
        FatalError("Matrix::relax: alpha must be positive");
    }

    const Count numCells = mesh_.numCells();

    if (phiPrevIter.size() != numCells)
    {
        FatalError
        (
            "Matrix::relax: phiPrevIter size mismatch "
            "with matrix size"
        );
    }

    // Store relaxation factor for setValues to recover
    // the pre-relaxation diagonal
    lastRelaxationFactor_ = alpha;

    const Scalar factor = (S(1.0) - alpha) / alpha;

    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar origDiag = cooValues_[diagOffset_ + cellIdx];

        // Scale diagonal: a_P <- a_P / alpha
        cooValues_[diagOffset_ + cellIdx] = origDiag / alpha;

        // Update RHS: b <- b + ((1-alpha)/alpha) * a_P * phiPrevIter
        vectorB_[cellIdx] += factor * origDiag * phiPrevIter[cellIdx];
    }
}


void Matrix::assemble()
{
    PETSC_CHECK(MatSetValuesCOO(matrixA_, cooValues_.data(), INSERT_VALUES));

    // Raw staging writes to vectorB_ bypass PETSc's object-state tracking,
    // so cached norms would go stale: VecNorm (ours and the KSP convergence
    // test's) would keep returning the first solve's |b|. The get/restore
    // pair is the public API for "the array changed" — restore bumps the
    // state; on wrapped storage it moves no data
    PetscScalar* rhsArray = nullptr;
    PETSC_CHECK(VecGetArray(rhsVec_, &rhsArray));
    PETSC_CHECK(VecRestoreArray(rhsVec_, &rhsArray));
}


void Matrix::explicitJacobiUpdate
(
    const ScalarField& phiOld,
    ScalarField& phiNew
) const
{
    const Count numCells = mesh_.numCells();

    // Single Jacobi sweep:
    //     phiNew_P = (b_P - sum_{N != P} A_PN phiOld_N) / A_PP.
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar diagonal = cooValues_[diagOffset_ + cellIdx];
        Scalar offDiagonalSum = S(0.0);

        for (const Index faceIdx : mesh_.cells()[cellIdx].faceIndices())
        {
            const Face& face = mesh_.faces()[faceIdx];

            if (face.isBoundary())
            {
                continue;
            }

            // Row cellIdx couples to the cell across each internal face
            const Index slot = faceSlot_[faceIdx];
            const bool isOwner = face.ownerCell() == cellIdx;
            const Index rowSlot = isOwner ? slot : slot + 1;
            const Index neighborIdx =
                isOwner ? face.neighborCell().value() : face.ownerCell();

            offDiagonalSum += cooValues_[rowSlot] * phiOld[neighborIdx];
        }

        phiNew[cellIdx] =
            (vectorB_[cellIdx] - offDiagonalSum) / (diagonal + vSmallValue);
    }
}

// ****************************** Private Methods *****************************

void Matrix::assembleInternalFace
(
    Index k,
    const Face& face,
    const TransportEquation& equation
)
{
    const Index ownerIdx = face.ownerCell();
    const Index neighborIdx = face.neighborCell().value();
    const Vector Sf = face.normal() * face.projectedArea();
    const Vector dPN =
        mesh_.cells()[neighborIdx].centroid()
      - mesh_.cells()[ownerIdx].centroid();
    const Scalar dPNMag = magnitude(dPN);
    const Vector ePN = dPN / (dPNMag + vSmallValue);

    // Orthogonal component (over-relaxed)
    const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePN)) * ePN;
    const Scalar Gammaf = equation.GammaFace[face.idx()];
    const Scalar aDiff = Gammaf * magnitude(Ef) / (dPNMag + vSmallValue);

    // Convection coefficients
    Scalar flowRate = S(0.0);
    Scalar aPConv = S(0.0);
    Scalar aNConv = S(0.0);

    if (equation.convection)
    {
        const ConvectionTerm& convection = *equation.convection;

        flowRate = convection.flowRate[face.idx()];

        // Implicit upwind split: a_P = max(F, 0), a_N = min(F, 0)
        aPConv = std::max(flowRate, S(0.0));
        aNConv = std::min(flowRate, S(0.0));
    }

    // Off-diagonal slots (exclusive to this face) and per-side diagonal
    // contributions gathered later by the owner/neighbor cells
    cooValues_[2 * k] = -aDiff + aNConv;            // A(owner, neighbor)
    cooValues_[2 * k + 1] = -aDiff - aPConv;        // A(neighbor, owner)
    faceDiag_[2 * k] = aDiff + aPConv;              // owner diagonal share
    faceDiag_[2 * k + 1] = aDiff - aNConv;          // neighbor diagonal share

    // Non-orthogonal correction (explicit)
    const Vector Tf = Sf - Ef;
    const Vector& gradPhiP = equation.gradPhi[ownerIdx];
    const Vector& gradPhiN = equation.gradPhi[neighborIdx];
    const Vector gradPhif =
        equation.gradScheme.faceGradient
        (
            equation.field,
            equation.phi,
            gradPhiP,
            gradPhiN,
            face.idx()
        );
    const Scalar nonOrthogonalFlux = Gammaf * dot(gradPhif, Tf);

    Scalar rhsOwner = nonOrthogonalFlux;
    Scalar rhsNeighbor = -nonOrthogonalFlux;

    // Deferred correction (explicit, convection only)
    if (equation.convection)
    {
        const ConvectionTerm& convection = *equation.convection;

        const Scalar deferredCorrection =
            convection.scheme.correction
            (
                face,
                equation.phi,
                gradPhiP,
                gradPhiN,
                flowRate
            );

        rhsOwner -= deferredCorrection;
        rhsNeighbor += deferredCorrection;
    }

    faceRhs_[2 * k] = rhsOwner;
    faceRhs_[2 * k + 1] = rhsNeighbor;
}


void Matrix::assembleBoundaryFace
(
    Index m,
    const Face& face,
    const TransportEquation& equation
)
{
    const Index ownerIdx = face.ownerCell();
    const BoundaryData& bc =
        bcManager_.fieldBC(face.patch()->get().patchName(), equation.field);
    const Vector Sf = face.normal() * face.projectedArea();
    const Vector ePf = normalized(face.dPf());
    const Scalar dPfMag = face.dPfMag();
    const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePf)) * ePf;
    const Scalar Gammaf = equation.GammaFace[face.idx()];
    const Scalar aDiff = Gammaf * magnitude(Ef) / (dPfMag + vSmallValue);
    const ConvectionTerm* convection =
        equation.convection ? &*equation.convection : nullptr;

    // Owner-cell contributions, gathered by the cell pass
    Scalar diag = S(0.0);
    Scalar rhs = S(0.0);

    using enum BCType;
    const BCType type = bc.type();

    switch (type)
    {
        case fixedValue:
        case noSlip:
        {
            // Dirichlet BC: phiB is prescribed
            Scalar phiB = S(0.0);

            if (type != noSlip)
            {
                phiB = bc.fixedScalarValue();
            }

            // Diffusion contribution
            diag += aDiff;
            rhs += aDiff * phiB;

            // Convection contribution
            if (convection != nullptr)
            {
                const Scalar aConv = convection->flowRate[face.idx()];

                rhs -= aConv * phiB;
            }

            break;
        }
        case fixedGradient:
        {
            const Scalar gradient = bc.fixedScalarGradient();
            const Scalar dn = dot(face.dPf(), face.normal());

            rhs += Gammaf * gradient * face.projectedArea();

            // Boundary value for convection: phi_b = phi_P + grad * dn
            if (convection != nullptr)
            {
                const Scalar aConv = convection->flowRate[face.idx()];

                diag += aConv;
                rhs -= aConv * gradient * dn;
            }

            break;
        }
        case zeroGradient:
        case kWallFunction:
        case nutWallFunction:
        case omegaWallFunction:
        {
            // Zero normal gradient: only convection
            if (convection != nullptr)
            {
                diag += convection->flowRate[face.idx()];
            }
            // No convection + zero gradient = no contribution
            break;
        }
        case symmetry:
        {
            // Symmetry plane carries zero normal mass flux
            if (!equation.velocity.has_value())
            {
                break;
            }

            const VelocityComponents& U = *equation.velocity;
            const Vector n = face.normal();
            const Scalar UxP = U.Ux[ownerIdx];
            const Scalar UyP = U.Uy[ownerIdx];
            const Scalar UzP = U.Uz[ownerIdx];

            Scalar ni = S(0.0);
            Scalar UnCross = S(0.0);    // sum_{j!=i} n_j * U_{P,j}

            switch (equation.field)
            {
                case Field::Ux:
                    ni = n.x();
                    UnCross = n.y() * UyP + n.z() * UzP;
                    break;
                case Field::Uy:
                    ni = n.y();
                    UnCross = n.x() * UxP + n.z() * UzP;
                    break;
                case Field::Uz:
                    ni = n.z();
                    UnCross = n.x() * UxP + n.y() * UyP;
                    break;
                default:
                    break;
            }

            diag += aDiff * ni * ni;
            rhs -= aDiff * ni * UnCross;
            break;
        }
        default:
        {
            // Unhandled BC type: default to zero gradient
            Warning
            (
                "Undefined boundary condition type for field "
              + Name(fieldToString(equation.field)) + " on patch "
              + face.patch()->get().patchName()
              + ". Applying zero gradient."
            );

            if (convection != nullptr)
            {
                diag += convection->flowRate[face.idx()];
            }
            break;
        }
    }

    boundaryDiag_[m] = diag;
    boundaryRhs_[m] = rhs;
}


void Matrix::setValues
(
    IndexListRef cellIndices,
    ScalarListRef values,
    ScalarListRef fractions
)
{
    const bool hasFractions = !fractions.empty();

    for (Index i = 0; i < cellIndices.size(); ++i)
    {
        const Index cellIdx = cellIndices[i];
        const Scalar f = hasFractions ? fractions[i] : S(1.0);

        if (f < rootSmallValue_)
        {
            continue;
        }

        const Scalar diag = cooValues_[diagOffset_ + cellIdx];

        if (f > S(1.0) - rootSmallValue_)
        {
            // Full constraint: zero the row and column couplings, moving
            // the known-value column coupling to the neighbors' RHS.
            // Every coupling of this cell lives on one of its internal
            // faces: slot = A(owner, neighbor), slot + 1 = the reverse.
            for (const Index faceIdx : mesh_.cells()[cellIdx].faceIndices())
            {
                const Face& face = mesh_.faces()[faceIdx];

                if (face.isBoundary())
                {
                    continue;
                }

                const Index slot = faceSlot_[faceIdx];
                const bool isOwner = face.ownerCell() == cellIdx;
                const Index rowSlot = isOwner ? slot : slot + 1;
                const Index colSlot = isOwner ? slot + 1 : slot;
                const Index neighborIdx =
                    isOwner ? face.neighborCell().value() : face.ownerCell();

                const Scalar coupling = cooValues_[colSlot];
                if (coupling != S(0.0))
                {
                    vectorB_[neighborIdx] -= coupling * values[i];
                    cooValues_[colSlot] = S(0.0);
                }

                cooValues_[rowSlot] = S(0.0);
            }

            vectorB_[cellIdx] = diag * values[i];
        }
        else
        {
            const Scalar diagPre =
                lastRelaxationFactor_ > S(0.0)
              ? lastRelaxationFactor_ * diag
              : diag;

            const Scalar coeff = f / (S(1.0) - f) * diagPre;

            cooValues_[diagOffset_ + cellIdx] += coeff;
            vectorB_[cellIdx] += coeff * values[i];
        }
    }
}
