# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
# SPDX-License-Identifier: Apache-2.0
#
# Catch2 + CTest unit and MPI test suite. Two executables:
#   turblyze_unit_tests - serial, never calls MPI_Init (asserts the
#                         uninitialized-MPI identity paths); per-TEST_CASE CTest
#                         entries via catch_discover_tests, labelled "serial".
#   turblyze_mpi_tests  - custom main (MPIRuntime + PETScRuntime), launched
#                         under mpirun at np=1/2/4 via plain add_test entries.
# Prerequisite: brew install catch2
# -----------------------------------------------------------------------------

find_package(Catch2 3 REQUIRED)
include(Catch)

# --- Shared fixture support (programmatic hex-box mesh, tolerance constants) ---
add_library(turblyze_test_support STATIC
    support/MeshFixtures.cpp
)

target_include_directories(turblyze_test_support PUBLIC support)

target_link_libraries(turblyze_test_support PUBLIC turblyze_core)

turblyze_apply_build_flags(turblyze_test_support)

# --- Serial unit tests (Catch2 stock main; no MPI, no PETSc ever) ---
add_executable(turblyze_unit_tests
    unit/Primitives/VectorTests.cpp
    unit/Primitives/TensorTests.cpp
    unit/Schemes/TimeSchemeTests.cpp
    unit/Schemes/InterpolationTests.cpp
    unit/Schemes/ConvectionSchemeTests.cpp
    unit/Schemes/GradientSchemeTests.cpp
    unit/Mesh/FaceGeometryTests.cpp
    unit/Mesh/CellGeometryTests.cpp
    unit/Mesh/MeshCountTests.cpp
    unit/Schemes/GradientLimiterTests.cpp
    unit/Fields/FieldTests.cpp
    unit/BoundaryConditions/BoundaryTypesTests.cpp
    unit/Models/TurbulenceStaticsTests.cpp
    unit/Case/CaseReaderTests.cpp
    unit/Case/RuntimeSelectionTests.cpp
    unit/Case/BoundaryFactoryTests.cpp
    unit/Parallel/CommTests.cpp
    unit/Parallel/ReduceSerialTests.cpp
    unit/Parallel/GlobalIndexTests.cpp
)

target_link_libraries(turblyze_unit_tests PRIVATE
    turblyze_test_support
    Catch2::Catch2WithMain
)

target_compile_definitions(turblyze_unit_tests PRIVATE
    TURBLYZE_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures"
)

turblyze_apply_build_flags(turblyze_unit_tests)

catch_discover_tests(turblyze_unit_tests
    PROPERTIES
        LABELS serial
        TIMEOUT 30
)

# --- MPI + PETSc tests (custom main in TestMain.cpp) ---
add_executable(turblyze_mpi_tests
    mpi/TestMain.cpp
    mpi/LinearSystem/MatrixRelaxTests.cpp
    mpi/LinearSystem/LinearSolverTests.cpp
    mpi/LinearSystem/MatrixConstraintTests.cpp
    mpi/Parallel/ReduceParallelTests.cpp
    mpi/Parallel/GlobalIndexParallelTests.cpp
    mpi/Parallel/HaloExchangeTests.cpp
    mpi/Parallel/ReduceConsistencyTests.cpp
    mpi/Models/TurbulenceBehaviorTests.cpp
)

target_link_libraries(turblyze_mpi_tests PRIVATE
    turblyze_test_support
    Catch2::Catch2
)

turblyze_apply_build_flags(turblyze_mpi_tests)

# One CTest entry per rank count; np=1 also goes through mpirun so prte
# behaviour is uniform. On macOS every launch carries the HWLOC_SYNTHETIC
# workaround (prte segfaults there without it); other platforms keep the
# real hardware topology.
foreach(np IN ITEMS 1 2 4)
    add_test(
        NAME turblyze_mpi_np${np}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                $<TARGET_FILE:turblyze_mpi_tests>
    )
    set_tests_properties(turblyze_mpi_np${np} PROPERTIES
        LABELS mpi
        PROCESSORS ${np}
        TIMEOUT 120
    )
    if(APPLE)
        set_property(TEST turblyze_mpi_np${np} PROPERTY
            ENVIRONMENT "HWLOC_SYNTHETIC=numa:1 core:8 pu:1"
        )
    endif()
endforeach()

# --- FatalError death probe (separate process; aborts on each trigger) ---
add_executable(turblyze_abort_probe probe/AbortProbe.cpp)

target_link_libraries(turblyze_abort_probe PRIVATE turblyze_core)

turblyze_apply_build_flags(turblyze_abort_probe)

# Each trigger must print "FATAL ERROR" and abort; the control must exit clean.
# The probe dies by std::abort (SIGABRT), which CTest would otherwise report as
# an abnormal termination, so the script owns the verdict: it asserts both the
# non-zero status and the message, and reports a clean exit as a failure.
foreach(trigger IN ITEMS unknownConvection normalizeZero)
    add_test(
        NAME turblyze_abort_${trigger}
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/probe/checkAbort.sh"
                $<TARGET_FILE:turblyze_abort_probe>
                ${trigger}
    )
    set_tests_properties(turblyze_abort_${trigger} PROPERTIES
        LABELS death
        TIMEOUT 30
    )
endforeach()

add_test(NAME turblyze_abort_control COMMAND turblyze_abort_probe ok)
set_tests_properties(turblyze_abort_control PROPERTIES
    LABELS death
    TIMEOUT 30
)

# A deterministic, fixed-iteration steady sphere run. The pinned Cd is
# machine-local: -march=native makes it reproducible only on the CPU it was
# measured on (default below: Apple M4). On another machine, re-pin from
# column 4 of the "Cd" line in regressionOut/reg_forces.txt via
# -DTURBLYZE_REGRESSION_CD_PIN=<value>, or widen the tolerance.
set(TURBLYZE_REGRESSION_CD_PIN "9.046208e-01" CACHE STRING
    "Pinned sphere Cd for the regression test (machine-local; re-pin per CPU)")
set(TURBLYZE_REGRESSION_CD_TOL "1e-4" CACHE STRING
    "Relative tolerance for the pinned sphere Cd")
set(TURBLYZE_SPHERE_MESH "${CMAKE_SOURCE_DIR}/inputFiles/sphere.msh")
if(NOT EXISTS "${TURBLYZE_SPHERE_MESH}")
    message(STATUS "Turblyze: sphere mesh absent - regression test disabled")
elseif(NOT TURBLYZE_USE_DOUBLE_PRECISION)
    message(STATUS "Turblyze: single precision - regression test disabled")
elseif(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Turblyze: not a Release build - regression test disabled")
elseif(NOT TURBLYZE_NATIVE_ARCH)
    message(STATUS "Turblyze: no -march=native - regression test disabled")
else()
    # Absolute mesh path, so any build directory works
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/regression/sphereRegression.in"
        "${CMAKE_CURRENT_BINARY_DIR}/sphereRegression"
        @ONLY
    )

    add_test(
        NAME turblyze_regression_sphere
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/regression/checkRegression.sh"
                $<TARGET_FILE:Turblyze>
                "${CMAKE_CURRENT_BINARY_DIR}/sphereRegression"
                "${CMAKE_BINARY_DIR}/regressionOut/reg_forces.txt"
                ${TURBLYZE_REGRESSION_CD_PIN}
                ${TURBLYZE_REGRESSION_CD_TOL}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
    set_tests_properties(turblyze_regression_sphere PROPERTIES
        LABELS regression
        TIMEOUT 300
    )
endif()
