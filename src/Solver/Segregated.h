/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Segregated.h
 * @brief Abstract base class for segregated pressure-correction algorithms
 *
 * @details Segregated owns the pressure-correction and Rhie-Chow
 * shared by every segregated coupling algorithm (SIMPLE, PISO):
 * implicit momentum assembly, the pressure-correction Poisson solve, the 
 * velocity and face-flux corrections, and the momentum diagonal coefficients.
 * Concrete algorithms derive from it and implement only
 * outerIteration() and algorithmName()
 *
 * @class Segregated
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <optional>

// Project headers
#include "MomentumTransport.h"
#include "ConvectionSchemes.h"
#include "Matrix.h"
#include "LinearSolvers.h"
#include "TransportEquation.h"

// ***************************** class Segregated *****************************

class Segregated : public MomentumTransport
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    Segregated
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const TimeScheme& timeScheme,
        const GradientScheme& gradScheme,
        const ConvectionSchemes& momentumConvectionScheme,
        LinearSolver& momentumSolver,
        LinearSolver& pressureSolver,
        TurbulenceModel& turbulence,
        const Vector& initialVelocity,
        Scalar initialPressure,
        Scalar deltaT,
        Scalar rho,
        Scalar mu,
        Scalar alphaU,
        Scalar alphaP,
        Count maxIterations,
        Scalar convergenceTolerance,
        Count nNonOrthogonalCorrectors,
        Count nOuterCorrectors,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    Segregated(const Segregated&) = delete;
    Segregated& operator=(const Segregated&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    Segregated(Segregated&&) = delete;
    Segregated& operator=(Segregated&&) = delete;

    /// Destructor
    ~Segregated() noexcept override = default;

// ***************************** Protected Members ****************************

protected:

    /// The convective face mass flux
    [[nodiscard]] const FaceFluxField& faceMassFlux() const noexcept override
    {
        return RhieChowFlowRate_;
    }

    /// Pressure correction RMS normalized by the pressure RMS
    [[nodiscard]] Scalar pressureResidual() const noexcept override;

    /// Build cell and face effective viscosity
    void updateEffectiveViscosity();

    /// Build the pressure-gradient source, velocity gradients and the
    /// transpose-gradient source; reset DU_ accumulators
    void assembleMomentum();

    /// Solve the three momentum components
    void solveMomentum(const TransientFields* prevStep);

    /// Compute the momentum diagonal coefficient DU_ = V / a_P
    void diagonalDU();

    /// Interpolate the face momentum diagonal DUf_ from DU_
    void buildFaceDiagonal();

    /// Update face mass fluxes using Rhie-Chow interpolation
    void updateRhieChowFlowRate(const TransientFields* prevStep);

    /// Assemble and solve the pressure correction equation
    void solvePressureCorrection();

    /// Apply velocity correction: U = U* - D*gradPCorr
    void correctVelocity();

    /// Update face mass fluxes from the pressure correction gradient
    void correctFlowRate();

    /// Update pressure with under-relaxation and reset pCorr
    void correctPressure();

    /// Add the transpose-gradient source to the momentum source terms
    void addTransposeGradientSource();

    /// Build the transient term for one velocity component
    [[nodiscard]] std::optional<TransientTerm> ddtTerm
    (
        const TransientFields* prevStep,
        ScalarField TransientFields::* phiPrevStep,
        ScalarField TransientFields::* ddtPrevStep
    ) const;

// **************************** Protected Accessors ***************************

    /// Momentum convection scheme
    [[nodiscard]] const ConvectionSchemes&
    momentumConvectionScheme() const noexcept
    {
        return momentumConvectionScheme_;
    }

    /// Matrix constructor and solver object
    [[nodiscard]] const Matrix& matrixConstruct() const noexcept
    {
        return matrixConstruct_;
    }
    [[nodiscard]] Matrix& matrixConstruct() noexcept
    {
        return matrixConstruct_;
    }

    /// Pressure gradient field
    [[nodiscard]] const VectorField& gradP() const noexcept
    {
        return gradP_;
    }
    [[nodiscard]] VectorField& gradP() noexcept
    {
        return gradP_;
    }

    /// Effective viscosity at face centres
    [[nodiscard]] const FaceData<Scalar>& nuEffFace() const noexcept
    {
        return nuEffFace_;
    }

    /// Momentum source terms
    [[nodiscard]] const ScalarField& UxSource() const noexcept
    {
        return UxSource_;
    }
    [[nodiscard]] const ScalarField& UySource() const noexcept
    {
        return UySource_;
    }
    [[nodiscard]] const ScalarField& UzSource() const noexcept
    {
        return UzSource_;
    }

// ****************************** Private Members *****************************

private:

// Segregated dependencies

    /// Reference to the momentum convection scheme
    const ConvectionSchemes& momentumConvectionScheme_;

    /// Linear solver for momentum equations
    LinearSolver& momentumSolver_;

    /// Linear solver for the pressure correction equation
    LinearSolver& pressureSolver_;

    /// Matrix constructor and solver object
    Matrix matrixConstruct_;

// Algorithm parameters

    /// Under-relaxation factor for velocity
    Scalar alphaU_;

    /// Under-relaxation factor for pressure
    Scalar alphaP_;

    /// Non-orthogonal corrector sub-iterations for the p' equation
    Count nNonOrthogonalCorrectors_;

// Pressure-correction fields

    /// Pressure correction field
    ScalarField pCorr_;

    /// Pressure correction gradient field
    VectorField gradPCorr_;

    /// Pressure gradient field
    VectorField gradP_;

// Face velocity fields

    /// Face velocity (current iteration)
    FaceData<Scalar> UxAvgf_;
    FaceData<Scalar> UyAvgf_;
    FaceData<Scalar> UzAvgf_;

    /// Face velocity (previous iteration)
    FaceData<Scalar> UxAvgPrevIterf_;
    FaceData<Scalar> UyAvgPrevIterf_;
    FaceData<Scalar> UzAvgPrevIterf_;

// Mass flux fields

    /// Mass flux through faces (Rhie-Chow)
    FaceFluxField RhieChowFlowRate_;

    /// Mass flux from the previous iteration
    FaceFluxField RhieChowFlowRatePrevIter_;

// Momentum diagonal coefficients

    /// Momentum diagonal coefficients
    ScalarField DU_;

    /// Face momentum diagonal coefficients
    FaceFluxField DUf_;

    /// Flag to allow one-time computation of DU_ per assembly
    bool DUComputed_ = false;

// Momentum assembly fields

    /// Effective viscosity (laminar + turbulent)
    ScalarField nuEff_;

    /// Effective viscosity at face centres
    FaceData<Scalar> nuEffFace_;

    /// Momentum source terms
    ScalarField UxSource_;
    ScalarField UySource_;
    ScalarField UzSource_;

// Pressure-correction assembly fields

    /// Mass imbalance source for the pressure correction equation
    ScalarField massImbalanceSrc_;

    /// Track pressure correction RMS before reset
    Scalar lastPressureCorrectionRMS_ = S(1e9);
};