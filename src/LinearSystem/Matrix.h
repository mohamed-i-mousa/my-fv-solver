/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Matrix.h
 * @brief Matrix assembly for finite volume discretization
 *
 * @details This header defines the Matrix class, which handles assembly of
 * sparse linear systems (Ax=b) for any transport equation type: momentum,
 * pressure correction, and turbulence, described by a TransportEquation
 * bundle (see TransportEquation.h).
 *
 * Coefficients are staged in a COO (coordinate) value array whose sparsity
 * pattern is fixed once at construction from the mesh topology and handed
 * to PETSc via MatSetPreallocationCOO. The slot layout is
 *
 *     [ 2 per internal face: off-diagonals | one per cell: diagonal ]
 *
 * Assembly needs three race-free passes:
 * Pass 1: face k touches only slots 2k/2k+1, no two faces collide.
 * Pass 2: face m touches only boundaryDiag_[m]/boundaryRhs_[m].
 * Pass 3: cell cellIdx reads shared scratch but writes only its own
 *         diagonal slot + RHS entry — one writer per cell.
 *
 * @class Matrix
 *
 * The Matrix class provides:
 * - Unified assembly for all transport equation types
 * - PETSc matrix/vector storage for the solve
 * - Deferred correction for higher-order convection schemes
 * - Non-orthogonal mesh corrections
 * - Implicit under-relaxation (Patankar) for solution stability
 * - Boundary condition integration during assembly
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <span>
#include <vector>

// External library headers
#include <petscmat.h>
#include <petscvec.h>

// Project headers
#include "Mesh.h"
#include "BoundaryConditions.h"
#include "Integer.h"
#include "TransportEquation.h"

// ******************************* class Matrix *******************************

class Matrix
{
public:

    // Type reductions for readability
    using ScalarListRef = std::span<const Scalar>;
    using ScalarSpan = std::span<Scalar>;

// ************************* Special Member Functions *************************

    /// Constructor: builds the PETSc system
    Matrix
    (
        const Mesh& mesh,
        const BoundaryConditions& boundaryConds
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    Matrix(Matrix&&) = delete;
    Matrix& operator=(Matrix&&) = delete;

    /// Destructor releases the PETSc matrix and RHS vector
    ~Matrix() noexcept;

// ****************************** Public Methods ******************************

    /// Build transport equation coefficients
    void buildMatrix(const TransportEquation& equation);

    /// Apply Patankar implicit under-relaxation to the staged values
    void relax(Scalar alpha, const ScalarField& phiPrevIter);

    /// Fix matrix rows to impose known cell values
    /// Replaces each constrained cell's equation with
    /// diag * phi[i] = diag * value[i] and moves the known-value coupling
    /// to unconstrained neighbors' RHS; call after relax(), before solve()
    void setValues
    (
        IndexListRef cellIndices,
        ScalarListRef values,
        ScalarListRef fractions = {}
    );

    /// Push the staged values into the PETSc matrix
    void assemble();

    /// One explicit Jacobi sweep of the staged system
    void explicitJacobiUpdate
    (
        const ScalarField& phiOld,
        ScalarField& phiNew
    ) const;

// ***************************** Accessor Methods *****************************

    /// Get the assembled PETSc matrix
    [[nodiscard]] Mat matrixA() const noexcept
    {
        return matrixA_;
    }

    /// Get the PETSc right-hand side vector
    [[nodiscard]] Vec rhsVec() const noexcept
    {
        return rhsVec_;
    }

    /// Get right-hand side vector b (const)
    [[nodiscard]] const ScalarList& vectorB() const noexcept
    {
        return vectorB_;
    }

    /// Get right-hand side vector b
    [[nodiscard]] ScalarList& vectorB() noexcept
    {
        return vectorB_;
    }

    /// Get the matrix diagonal (const)
    [[nodiscard]] ScalarListRef diagonal() const noexcept
    {
        return {cooValues_.data() + diagOffset_, mesh_.numCells()};
    }

    /// Get the staged matrix diagonal
    [[nodiscard]] ScalarSpan diagonal() noexcept
    {
        return {cooValues_.data() + diagOffset_, mesh_.numCells()};
    }

// ****************************** Private Members *****************************

private:

    /// Mesh data reference
    const Mesh& mesh_;

    /// Boundary data references
    const BoundaryConditions& bcManager_;

    /// PETSc sparse matrix
    Mat matrixA_ = nullptr;

    /// PETSc RHS vectorB_ storage
    Vec rhsVec_ = nullptr;

    /// Staged COO values: [internal-face off-diagonals | cell diagonals]
    ScalarList cooValues_;

    /// Right-hand side staging
    ScalarList vectorB_;

    /// Per-face diagonal contributions [2 per internal face: owner, neighbor]
    ScalarList faceDiag_;

    /// Per-face RHS contributions [2 per internal face: owner, neighbor]
    ScalarList faceRhs_;

    /// Per-boundary-face diagonal contribution to the owner cell
    ScalarList boundaryDiag_;

    /// Per-boundary-face RHS contribution to the owner cell
    ScalarList boundaryRhs_;

    /// Face indices of every internal face, in slot order
    IndexList internalFaces_;

    /// Face indices of every boundary face, in slot order
    IndexList boundaryFaces_;

    /// faceIdx -> first off-diagonal slot (internal) or boundary ordinal
    IndexList faceSlot_;

    /// First diagonal slot in cooValues_ (= 2 * number of internal faces)
    Index diagOffset_ = 0;

    /// Relaxation factor from last relax() call (0 = not relaxed)
    Scalar lastRelaxationFactor_ = S(0.0);

    /// Threshold below which f/(1-f) overflows
    inline static const Scalar rootSmallValue_ = std::sqrt(smallValue);


// ****************************** Private Methods *****************************

private:

    /// Stage internal-face coefficients into slot pair 2k/2k+1
    void assembleInternalFace
    (
        Index k,
        const Face& face,
        const TransportEquation& equation
    );

    /// Stage boundary-face contributions into boundary scratch slot m
    void assembleBoundaryFace
    (
        Index m,
        const Face& face,
        const TransportEquation& equation
    );

};
