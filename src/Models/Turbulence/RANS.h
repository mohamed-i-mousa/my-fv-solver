/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file RANS.h
 * @brief Abstract base for two-equation Reynolds-Averaged Navier-Stokes models
 *
 * @details RANS sits between the root TurbulenceModel interface and concrete
 * two-equation models. It owns the state and services common to such models:
 * mesh/BC views, equation assembly services, turbulent viscosity, k transport
 * state, wall-function geometry/diagnostics, and residual bookkeeping.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <optional>
#include <utility>
#include <vector>

// Project headers
#include "Scalar.h"
#include "Integer.h"
#include "StringTypes.h"
#include "Vector.h"
#include "CellData.h"
#include "FaceData.h"
#include "TransportEquation.h"
#include "TurbulenceModel.h"

// *************************** Forward Declarations ***************************

class Mesh;
class BoundaryConditions;
class TimeScheme;
class GradientScheme;
class ConvectionScheme;
class LinearSolver;
class Matrix;
class Face;
enum class BCType;
enum class Field;

// ******************************** class RANS ********************************

class RANS : public TurbulenceModel
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    RANS
    (
        const Mesh& mesh,
        const BoundaryConditions& bc,
        const TimeScheme& timeScheme,
        const GradientScheme& gradScheme,
        const ConvectionScheme& kScheme,
        LinearSolver& kSolver,
        const ConvectionScheme& dissipationScheme,
        LinearSolver& dissipationSolver,
        Scalar deltaT,
        Scalar nu,
        Scalar alphaK,
        Scalar alphaDissipation,
        bool debug
    );

    /// Copy constructor and assignment - Not copyable (reference members)
    RANS(const RANS&) = delete;
    RANS& operator=(const RANS&) = delete;

    /// Move constructor and assignment - Not movable (reference members)
    RANS(RANS&&) = delete;
    RANS& operator=(RANS&&) = delete;

    /// Destructor - out-of-line to anchor the unique_ptr<Matrix> member
    ~RANS() noexcept override;

// ***************************** Turbulence Solve *****************************

    /// Solve turbulence equations for current iteration
    void solve
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz,
        const FaceFluxField& flowRateFace,
        const TensorField& gradU
    ) override = 0;

    /// k and dissipation at the start of a transient time step
    void beginTimeStep() override;

    /// Roll the Crank-Nicolson stored time derivatives forward one step
    void updatePrevStepDerivatives() override;

// ************************ Inlet Condition Calculators ***********************

    /// Calculate inlet/initial turbulent kinetic energy
    [[nodiscard]] static Scalar inletK
    (
        const Vector& velocity,
        Scalar turbulenceIntensity
    ) noexcept;

// ***************************** Accessor Methods *****************************

    /// Get turbulent kinematic viscosity field
    [[nodiscard]] const ScalarField&
    turbulentViscosity() const noexcept override
    {
        return nut_;
    }

    /// Get turbulent viscosity for a boundary face
    [[nodiscard]] Scalar boundaryTurbulentViscosity
    (
        const Face& face,
        const BoundaryConditions& bcManager
    ) const override;

    /// Wall shear stress (tau/rho) from the wall-function state and velocity
    [[nodiscard]] FaceData<Scalar> wallShearStress
    (
        const ScalarField& Ux,
        const ScalarField& Uy,
        const ScalarField& Uz
    ) const override;

    /// Get turbulent kinetic energy field
    [[nodiscard]] const ScalarField& k() const noexcept
    {
        return k_;
    }

    /// Get the model's dissipation field
    [[nodiscard]] virtual const ScalarField& dissipation() const noexcept = 0;

    /// Name of the dissipation field for output labelling ("omega"/"epsilon")
    [[nodiscard]] virtual Name dissipationName() const noexcept = 0;

    /// Whether wall-distance initialization converged
    [[nodiscard]] bool wallDistanceConverged() const noexcept override
    {
        return wallDistanceConverged_;
    }

    /// Cell-centered scalar fields exported by RANS models (name, field)
    [[nodiscard]] CellDataPair cellDataOutputs() const override;

    /// Boundary scalar fields exported by RANS models (name, field)
    [[nodiscard]] BoundaryDataPair boundaryDataOutputs() const override;

    /// Residuals contributed by RANS models to convergence checks (name, value)
    [[nodiscard]] ResidualPair residualOutputs() const override;


protected:

// **************************** Protected Accessors ***************************

// Dependencies and services

    /// Mesh view (nodes, faces, cells)
    [[nodiscard]] const Mesh& mesh() const noexcept { return mesh_; }

    /// Boundary conditions view
    [[nodiscard]] const BoundaryConditions& bcManager() const noexcept
    {
        return bcManager_;
    }

    /// Gradient scheme
    [[nodiscard]] const GradientScheme& gradientScheme() const noexcept
    {
        return gradientScheme_;
    }

    /// Matrix constructor (mutable Matrix through a const unique_ptr)
    [[nodiscard]] const std::unique_ptr<Matrix>&
    matrixConstruct() const noexcept
    {
        return matrixConstruct_;
    }

    /// k convection scheme
    [[nodiscard]] const ConvectionScheme& kConvectionScheme() const noexcept
    {
        return kConvectionScheme_;
    }

    /// Linear solver for the k equation
    [[nodiscard]] const LinearSolver& kSolver() const noexcept
    {
        return kSolver_;
    }
    [[nodiscard]] LinearSolver& kSolver() noexcept { return kSolver_; }

    /// Dissipation convection scheme
    [[nodiscard]] const ConvectionScheme&
    dissipationConvectionScheme() const noexcept
    {
        return dissipationConvectionScheme_;
    }

    /// Linear solver for the dissipation equation
    [[nodiscard]] const LinearSolver& dissipationSolver() const noexcept
    {
        return dissipationSolver_;
    }
    [[nodiscard]] LinearSolver& dissipationSolver() noexcept
    {
        return dissipationSolver_;
    }

// Physical and algorithm parameters

    /// Laminar kinematic viscosity
    [[nodiscard]] Scalar nu() const noexcept { return nu_; }

    /// Under-relaxation factor for the k equation
    [[nodiscard]] Scalar alphaK() const noexcept { return alphaK_; }

    /// Under-relaxation factor for the dissipation equation
    [[nodiscard]] Scalar alphaDissipation() const noexcept
    {
        return alphaDissipation_;
    }

    /// Whether verbose console output is enabled
    [[nodiscard]] bool debug() const noexcept { return debug_; }

// Common transport state

    /// Turbulent kinematic viscosity (native short-name accessor)
    [[nodiscard]] const ScalarField& nut() const noexcept { return nut_; }
    [[nodiscard]] ScalarField& nut() noexcept { return nut_; }

    /// Mutable turbulent kinetic energy (the const overload is public)
    [[nodiscard]] ScalarField& k() noexcept { return k_; }

    /// Previous-iteration k snapshot
    [[nodiscard]] const ScalarField& kPrev() const noexcept { return kPrev_; }
    [[nodiscard]] ScalarField& kPrev() noexcept { return kPrev_; }

    /// k from the previous time step (phi^n)
    [[nodiscard]] const ScalarField& kPrevStep() const noexcept
    {
        return kPrevStep_;
    }

    /// Dissipation from the previous time step (phi^n)
    [[nodiscard]] const ScalarField& dissipationPrevStep() const noexcept
    {
        return dissipationPrevStep_;
    }

    /// Stored old time derivative of k for Crank-Nicolson
    [[nodiscard]] const ScalarField& kDdtPrevStep() const noexcept
    {
        return kDdtPrevStep_;
    }

    /// Stored old time derivative of dissipation for Crank-Nicolson
    [[nodiscard]] const ScalarField& dissipationDdtPrevStep() const noexcept
    {
        return dissipationDdtPrevStep_;
    }

    /// Cell gradient of k
    [[nodiscard]] const VectorField& gradK() const noexcept { return gradK_; }
    [[nodiscard]] VectorField& gradK() noexcept { return gradK_; }

// Wall distance and wall-function state

    /// Distance to nearest wall
    [[nodiscard]] const ScalarField& wallDistance() const noexcept
    {
        return wallDistance_;
    }

    /// Owner-cell to wall-face perpendicular distance
    [[nodiscard]] const FaceData<Scalar>& y() const noexcept { return y_; }

    /// y+ on wall-function faces
    [[nodiscard]] const FaceData<Scalar>& yPlus() const noexcept
    {
        return yPlus_;
    }

    /// Wall-function nut values on wall faces
    [[nodiscard]] const FaceData<Scalar>& nutWall() const noexcept
    {
        return nutWall_;
    }
    [[nodiscard]] FaceData<Scalar>& nutWall() noexcept { return nutWall_; }

    /// Area-based weight per wall face
    [[nodiscard]] const FaceData<Scalar>& wallFaceWeight() const noexcept
    {
        return wallFaceWeight_;
    }

    /// Indices of faces with model wall-function BCs
    [[nodiscard]] const IndexList& wallFunctionFaceIndices() const noexcept
    {
        return wallFunctionFaceIndices_;
    }

    /// Unique cell indices adjacent to wall-function faces
    [[nodiscard]] const IndexList& wallCellIndices() const noexcept
    {
        return wallCellIndices_;
    }

    /// Wall-to-total boundary area fraction per wall cell
    [[nodiscard]] const ScalarList& wallCellFraction() const noexcept
    {
        return wallCellFraction_;
    }

    /// y+ crossover between viscous sublayer and log region
    [[nodiscard]] Scalar yPlusLam() const noexcept { return yPlusLam_; }

// ****************************** Shared Methods ******************************

    /// Build the transient term for one field, or nullopt if steady
    [[nodiscard]] std::optional<TransientTerm> ddtTerm
    (
        const ScalarField& phiPrevStep,
        const ScalarField& ddtPrevStep
    ) const;

    /// Map cell-centered diffusion coefficients to faces for assembly
    void cellToFaceDiffusion
    (
        const ScalarField& cellGamma,
        FaceFluxField& faceGamma
    ) const;

    /// Update wall distance field using mesh-wave coordinate propagation
    void updateWallDistance();

    /// Compute y+ crossover via fixed-point iteration
    void updateYPlusLam(Scalar kappa, Scalar E);

    /// Build wall-function face lists and area weights
    void initializeWallFunctionGeometry
    (
        const BoundaryConditions& bcManager,
        Field wallFunctionField,
        BCType wallFunctionType
    );

    /// Model-specific Cμ^0.25 used by the wall functions
    [[nodiscard]] virtual Scalar cmu25() const noexcept = 0;

    /// Update y+ field on wall-function faces
    void updateYPlus();

    /// Compute strain-rate magnitude: ||S|| = sqrt(2 S_ij S_ij)
    [[nodiscard]] ScalarField computeStrainRateMagnitude
    (
        const TensorField& gradU
    ) const;

    /// Compute cell velocity divergence, allocated on demand
    [[nodiscard]] ScalarField velocityDivergence
    (
        const FaceFluxField& flowRateFace
    ) const;

    /// Update both residuals against the pre-solve snapshots
    void updateResiduals
    (
        const ScalarField& dissipationField,
        const ScalarField& dissipationPrevField
    );

    /// Compute normalised field change against a previous snapshot
    [[nodiscard]] Scalar normalisedFieldResidual
    (
        const ScalarField& field,
        const ScalarField& previousField
    ) const;

// ****************************** Private Members *****************************

private:

// Dependencies and services

    /// Mesh view (nodes, faces, cells)
    const Mesh& mesh_;

    /// Reference to BCs
    const BoundaryConditions& bcManager_;

    /// Time-derivative discretization scheme
    const TimeScheme& timeScheme_;

    /// Reference to gradient scheme
    const GradientScheme& gradientScheme_;

    /// Matrix constructor
    std::unique_ptr<Matrix> matrixConstruct_;

    /// Reference to k convection scheme
    const ConvectionScheme& kConvectionScheme_;

    /// Linear solver for k equation
    LinearSolver& kSolver_;

    /// Reference to dissipation convection scheme
    const ConvectionScheme& dissipationConvectionScheme_;

    /// Linear solver for dissipation equation
    LinearSolver& dissipationSolver_;

// Physical and algorithm parameters

    /// Laminar kinematic viscosity
    Scalar nu_;

    /// Time step size [s] (transient runs)
    Scalar deltaT_;

    /// Under-relaxation factor for k equation
    Scalar alphaK_;

    /// Under-relaxation factor for dissipation equation
    Scalar alphaDissipation_;

    /// Enable verbose console output
    bool debug_;

// Common transport state

    /// Turbulent kinematic viscosity
    ScalarField nut_{S(0.0)};

    /// Turbulent kinetic energy
    ScalarField k_{S(1e-6)};

    /// Previous-iteration k snapshot for residual computation
    ScalarField kPrev_;

    /// k from the previous time step (phi^n) for the transient term
    ScalarField kPrevStep_;

    /// Dissipation from the previous time step (phi^n)
    ScalarField dissipationPrevStep_;

    /// Stored old time derivative of k for Crank-Nicolson
    ScalarField kDdtPrevStep_;

    /// Stored old time derivative of dissipation for Crank-Nicolson
    ScalarField dissipationDdtPrevStep_;

    /// Cell gradient of k
    VectorField gradK_;

    /// Normalised k change from the most recent solve
    Scalar lastKResidual_ = S(1e9);

    /// Normalised dissipation change from the most recent solve
    Scalar lastDissipationResidual_ = S(1e9);

// Wall distance and wall-function state

    /// Distance to nearest wall
    ScalarField wallDistance_{S(1.0)};

    /// Coordinates of the nearest wall point per cell (for mesh-wave)
    VectorField nearestWallPoint_;

    /// meshWave wall-distance loop convergence flag
    bool wallDistanceConverged_ = false;

    /// Owner-cell to wall-face perpendicular distance
    FaceData<Scalar> y_;

    /// y+
    FaceData<Scalar> yPlus_;

    /// Wall-function nut values on wall faces
    FaceData<Scalar> nutWall_;

    /// Area-based weight per wall face (face area / total wall area of cell)
    FaceData<Scalar> wallFaceWeight_;

    /// Indices into mesh_.faces for faces with model wall-function BCs
    IndexList wallFunctionFaceIndices_;

    /// Unique cell indices adjacent to wall-function faces
    IndexList wallCellIndices_;

    /// Wall-to-total boundary area fraction per wall cell
    ScalarList wallCellFraction_;

    /// y+ crossover between viscous sublayer and log region
    Scalar yPlusLam_ = S(11.225);
};
