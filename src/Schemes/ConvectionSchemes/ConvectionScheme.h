/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ConvectionScheme.h
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

// ************************** class ConvectionScheme **************************

class ConvectionScheme
{
public:

// ************************* Special Member Functions *************************

    /// Copy constructor and assignment - Slicing problem
    ConvectionScheme(const ConvectionScheme&) = delete;
    ConvectionScheme& operator=(const ConvectionScheme&) = delete;

    /// Move constructor and assignment - Slicing problem
    ConvectionScheme(ConvectionScheme&&) = delete;
    ConvectionScheme& operator=(ConvectionScheme&&) = delete;

    /// Default destructor
    virtual ~ConvectionScheme() = default;

// **************************** Runtime Selection ****************************

    /// Construct the convection scheme selected by name
    [[nodiscard]] static std::unique_ptr<ConvectionScheme> create
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
    ConvectionScheme() = default;
};
