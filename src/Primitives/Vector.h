/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Vector.h
 * @brief 3D vector class for geometric and mathematical operations in CFD
 *
 * @details This header defines a 3D vector class that serves as the foundation
 * for all vector-based calculations in the CFD solver. The Vector class
 * provides essential mathematical operations required in the finite volume
 * discretization and mesh operations.
 *
 * All member and free functions are HD-annotated so the same inline math is
 * callable from CUDA kernels; error handling degrades to assert on the
 * device pass (see HostDevice.h). Stream output stays host-only.
 *
 * @class Vector
 * - Components access and manipulation (x, y, z coordinates)
 * - Arithmetic operations (addition, subtraction, scalar multiplication)
 * - Vector operations (dot product, cross product, normalization)
 * - Equality comparison operator
 * - Stream I/O operators for debugging
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <cassert>
#include <cmath>
#include <iosfwd>

// Project headers
#include "HostDevice.h"
#include "Scalar.h"
#include "ErrorHandler.h"

// ******************************* class Vector *******************************

class Vector
{
public:

/// ************************* Special Member Functions ************************

    /// Default constructor
    Vector() noexcept = default;

    /// Construct vector with specified components
    HD Vector(Scalar xValue, Scalar yValue, Scalar zValue) noexcept
    :
        x_(xValue),
        y_(yValue),
        z_(zValue)
    {}

// ****************************** Setter Methods ******************************

    /// Set X component
    HD void setX(Scalar xValue) noexcept { x_ = xValue; }

    /// Set Y component
    HD void setY(Scalar yValue) noexcept { y_ = yValue; }

    /// Set Z component
    HD void setZ(Scalar zValue) noexcept { z_ = zValue; }

// ***************************** Accessor Methods *****************************

    /// Get X component
    [[nodiscard]] HD Scalar x() const noexcept { return x_; }

    /// Get Y component
    [[nodiscard]] HD Scalar y() const noexcept { return y_; }

    /// Get Z component
    [[nodiscard]] HD Scalar z() const noexcept { return z_; }

// ***************************** Operator Methods *****************************

    /// Vector addition operator
    HD Vector operator+(const Vector& other) const noexcept
    {
        Vector result(*this);
        result += other;
        return result;
    }

    /// Vector subtraction operator
    HD Vector operator-(const Vector& other) const noexcept
    {
        Vector result(*this);
        result -= other;
        return result;
    }

    /// Scalar multiplication operator
    HD Vector operator*(Scalar scalar) const noexcept
    {
        Vector result(*this);
        result *= scalar;
        return result;
    }

    /// Scalar division operator
    HD Vector operator/(Scalar scalar) const noexcept
    {
        Vector result(*this);
        result /= scalar;
        return result;
    }

    /// Compound addition assignment operator
    HD Vector& operator+=(const Vector& other) noexcept
    {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;

        return *this;
    }

    /// Compound subtraction assignment operator
    HD Vector& operator-=(const Vector& other) noexcept
    {
        x_ -= other.x_;
        y_ -= other.y_;
        z_ -= other.z_;

        return *this;
    }

    /// Compound multiplication assignment operator
    HD Vector& operator*=(Scalar scalar) noexcept
    {
        x_ *= scalar;
        y_ *= scalar;
        z_ *= scalar;

        return *this;
    }

    /// Compound division assignment operator
    HD Vector& operator/=(Scalar scalar) noexcept
    {
#ifdef __CUDA_ARCH__
        // Device code cannot throw: the guard degrades to an assert
        assert(std::abs(scalar) > vSmallValue);
#else
        if (std::abs(scalar) <= vSmallValue)
        {
            FatalError("Division by zero in Vector::operator/=");
        }
#endif

        const Scalar inverse = S(1.0) / scalar;
        x_ *= inverse;
        y_ *= inverse;
        z_ *= inverse;

        return *this;
    }

    /// Equality comparison operator
    HD bool operator==(const Vector& other) const noexcept
    {
        return (std::abs(x_ - other.x_) <= vSmallValue)
            && (std::abs(y_ - other.y_) <= vSmallValue)
            && (std::abs(z_ - other.z_) <= vSmallValue);
    }

// ****************************** Private Members *****************************

private:

    /// x, y, z components of the vector
    Scalar x_ = S(0.0);
    Scalar y_ = S(0.0);
    Scalar z_ = S(0.0);
};

// *************************** Non-Member Functions ***************************

/// Compute dot product of two vectors
[[nodiscard]] HD inline Scalar dot(const Vector& p1, const Vector& p2) noexcept
{
    return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}

/// Compute cross product of two vectors
[[nodiscard]] HD inline Vector cross
(
    const Vector& p1,
    const Vector& p2
) noexcept
{
    return
        Vector
        (
            p1.y() * p2.z() - p1.z() * p2.y(),
            p1.z() * p2.x() - p1.x() * p2.z(),
            p1.x() * p2.y() - p1.y() * p2.x()
        );
}

/// Squared magnitude of a vector
[[nodiscard]] HD inline Scalar magnitudeSquared(const Vector& v) noexcept
{
    return v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
}

/// Magnitude of a vector
[[nodiscard]] HD inline Scalar magnitude(const Vector& v) noexcept
{
    return std::sqrt(magnitudeSquared(v));
}

/// Return a normalized copy of a vector
[[nodiscard]] HD inline Vector normalized(const Vector& v) noexcept
{
    const Scalar mag = magnitude(v);

#ifdef __CUDA_ARCH__
    // Device code cannot throw: the guard degrades to an assert
    assert(mag >= vSmallValue);
#else
    if (mag < vSmallValue)
    {
        FatalError("Division by zero in normalized(Vector)");
    }
#endif

    const Scalar inverse = S(1.0) / mag;
    return Vector(v.x() * inverse, v.y() * inverse, v.z() * inverse);
}

/// Scalar multiplication operator
HD inline Vector operator*(Scalar scalar, const Vector& p) noexcept
{
    return p * scalar;
}

/// Stream output operator for Vector
std::ostream& operator<<(std::ostream& os, const Vector& p);
