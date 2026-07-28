/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file AbortProbe.cpp
 * @brief Subprocess that deterministically triggers a FatalError
 *
 * @details FatalError is [[noreturn]] and terminates via std::abort (MPI is
 * never initialised here), so it cannot be caught in-process. This tiny binary
 * hits one FatalError path chosen by argv[1]; the CTest harness runs it as a
 * separate process and matches "FATAL ERROR" in its output. The "ok" mode is
 * a control that returns cleanly, proving the probe does not always abort.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <string_view>

// Project headers
#include "ConvectionScheme.h"
#include "Vector.h"
#include "Scalar.h"

// ********************************* Entry Point ******************************

int main(int argc, char** argv)
{
    const std::string_view mode = argc > 1 ? argv[1] : std::string_view{};

    // Control path: constructs nothing that aborts
    if (mode == "ok")
    {
        return 0;
    }

    // Unknown runtime selection funnels through FatalError
    if (mode == "unknownConvection")
    {
        const auto scheme = ConvectionScheme::create("NotAScheme");
        (void)scheme;   // unreachable: create aborts first
        return 0;
    }

    // Normalising a zero vector divides by zero and aborts
    if (mode == "normalizeZero")
    {
        const Vector unit = normalized(Vector(S(0.0), S(0.0), S(0.0)));
        (void)unit;     // unreachable: normalized aborts first
        return 0;
    }

    // Unknown mode: neither a clean exit nor a FatalError
    return 2;
}