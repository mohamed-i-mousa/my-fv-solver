/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TransportEquation.h
 * @brief Data bundle describing a scalar transport equation
 *
 * @details This header defines TransientTerm, ConvectionTerm,
 * VelocityComponents and TransportEquation. TransportEquation bundles all data
 * describing a scalar transport equation (field, transient, convection,
 * diffusion, source, velocity, gradients) needed by Matrix::buildMatrix().
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <optional>

// Project headers
#include "CellData.h"
#include "Field.h"
#include "FaceData.h"
#include "ConvectionScheme.h"
#include "GradientScheme.h"

// *************************** Forward Declarations ***************************

class TimeScheme;

// *************************** struct TransientTerm ***************************

struct TransientTerm
{
    /// Time-derivative discretization scheme
    const TimeScheme& scheme;

    /// Time step size
    Scalar deltaT;

    /// Field values at the previous time step (phi^n)
    const ScalarField& phiPrevStep;

    /// Stored old time derivative for Crank-Nicolson (nullptr if unused)
    const ScalarField* ddtPrevStep;
};

// *************************** struct ConvectionTerm **************************

struct ConvectionTerm
{
    /// Face volumetric flow rates
    const FaceFluxField& flowRate;

    /// Convection discretization scheme
    const ConvectionScheme& scheme;
};

// ************************* struct VelocityComponents ************************

struct VelocityComponents
{
    /// Velocity x-component field
    const ScalarField& Ux;

    /// Velocity y-component field
    const ScalarField& Uy;

    /// Velocity z-component field
    const ScalarField& Uz;
};

// ************************* struct TransportEquation *************************

struct TransportEquation
{

// *********************************** Field **********************************

    /// Field identifier
    Field field;

    /// Current cell-centered field values (mutable for zero-copy solve)
    ScalarField& phi;

// ******************************** Transient *********************************

    /// Complete transient term d(phi)/dt discretization (nullopt = steady)
    std::optional<TransientTerm> transient = std::nullopt;

// ******************************** Convection ********************************

    /// Complete convection term div(F * phi) (nullopt = no convection)
    std::optional<ConvectionTerm> convection = std::nullopt;

// ********************************* Diffusion ********************************

    /// Pre-interpolated face diffusion coefficient for div(Gamma * grad(phi))
    const FaceFluxField& GammaFace;

// ********************************** Source **********************************

    /// Explicit source term field
    const ScalarField& source;

// **************************** Velocity Components ***************************

    /// Velocity components for vector boundary conditions (symmetry).
    std::optional<VelocityComponents> velocity = std::nullopt;

// ************************** Gradient Reconstruction *************************

    /// Pre-computed cell gradients of phi
    const VectorField& gradPhi;

    /// Gradient reconstruction scheme
    const GradientScheme& gradScheme;
};
