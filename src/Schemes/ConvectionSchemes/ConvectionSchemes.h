/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ConvectionSchemes.h
 * @brief Abstract base class for convection discretization schemes
 *
 * @details This header defines the polymorphic convection scheme interface for
 * discretizing convective fluxes in transport equations. All schemes use
 * stable first-order upwind coefficients in the implicit matrix. Higher-order
 * schemes add an explicit deferred correction term through correction().
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <memory>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"
#include "Vector.h"
#include "Face.h"
#include "CellData.h"

// ************************** class ConvectionSchemes *************************

class ConvectionSchemes
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment - Slicing problem
    ConvectionSchemes(const ConvectionSchemes&) = delete;
    ConvectionSchemes& operator=(const ConvectionSchemes&) = delete;

    /// Move constructor and assignment - Slicing problem
    ConvectionSchemes(ConvectionSchemes&&) = delete;
    ConvectionSchemes& operator=(ConvectionSchemes&&) = delete;

    /// Default destructor
    virtual ~ConvectionSchemes() = default;

// **************************** Runtime Selection ****************************

    /// Construct the convection scheme selected by name
    [[nodiscard]] static std::unique_ptr<ConvectionSchemes> create
    (
        const Name& schemeName
    );

    /// Names of every selectable convection scheme
    [[nodiscard]] static NameList availableSchemes();

// ****************************** Public Methods ******************************

    /// Calculate higher-order deferred correction term
    [[nodiscard]] virtual Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const = 0;

// ***************************** Protected Methods ****************************

protected:

    /// Default constructor
    ConvectionSchemes() = default;
};
