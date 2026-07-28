# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)
project(Turblyze VERSION "1.0.0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Build options
option(TURBLYZE_USE_DOUBLE_PRECISION
       "Build with double-precision Scalar (otherwise float)" ON)
option(TURBLYZE_NATIVE_ARCH
       "Optimize for the build host CPU (-march=native)" ON)

# Find Eigen3 package (for least-squares precompute only)
find_package(Eigen3 REQUIRED NO_MODULE)

# Find PETSc via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

# Find MPI (initialized through PETSc; called directly in src/Parallel/)
find_package(MPI REQUIRED COMPONENTS CXX)

# Find METIS (mesh partitioning, used by src/Parallel/MeshDecomposer)
find_path(METIS_INCLUDE_DIR metis.h REQUIRED)
find_library(METIS_LIBRARY metis REQUIRED)

# Find the HDF5 C library
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED COMPONENTS C)

# Warning flags, shared by the solver library and the test targets: warning-
# clean under this set is part of the Definition of Done (GCC/Clang/AppleClang)
set(TURBLYZE_WARNING_FLAGS
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wsign-conversion
    -Wdouble-promotion
    -Wformat=2
    -Wcast-align
    -Wunused
    -Wnull-dereference
    -Wimplicit-fallthrough
    -Wold-style-cast
    -Woverloaded-virtual
    -Wnon-virtual-dtor
)

# Release optimization flags, likewise shared with the test targets
set(TURBLYZE_RELEASE_FLAGS -O3)
if(TURBLYZE_NATIVE_ARCH)
    list(APPEND TURBLYZE_RELEASE_FLAGS -march=native)
endif()

# One definition of the shared warning and Release flags for every target
# (solver core, executable, and the test targets in tests/)
function(turblyze_apply_build_flags target)
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:${TURBLYZE_WARNING_FLAGS}>
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:${TURBLYZE_RELEASE_FLAGS}>
        )
    endif()
endfunction()

# --- Solver core: every source except main.cpp, so the test targets can link
#     the exact same objects the executable does (OBJECT library keeps the
#     np=1 behaviour bit-identical to the pre-split single-executable build) ---
add_library(turblyze_core OBJECT
    src/Primitives/Vector.cpp
    src/Primitives/Tensor.cpp
    src/Primitives/Logger.cpp
    src/Mesh/Face.cpp
    src/Mesh/Cell.cpp
    src/Mesh/MeshReader.cpp
    src/Mesh/MeshCreator.cpp
    src/Mesh/MeshChecker.cpp
    src/Fields/Field.cpp
    src/BoundaryConditions/BoundaryConditions.cpp
    src/BoundaryConditions/BCLoader.cpp
    src/BoundaryConditions/BoundaryTypes/BoundaryType.cpp
    src/BoundaryConditions/BoundaryTypes/FixedValue.cpp
    src/BoundaryConditions/BoundaryTypes/NoSlip.cpp
    src/BoundaryConditions/BoundaryTypes/FixedGradient.cpp
    src/BoundaryConditions/BoundaryTypes/ZeroGradient.cpp
    src/BoundaryConditions/BoundaryTypes/WallFunction.cpp
    src/BoundaryConditions/BoundaryTypes/Symmetry.cpp
    src/Schemes/GradientSchemes/GradientScheme.cpp
    src/Schemes/GradientSchemes/LeastSquares.cpp
    src/Schemes/ConvectionSchemes/ConvectionScheme.cpp
    src/Schemes/ConvectionSchemes/UpwindScheme.cpp
    src/Schemes/ConvectionSchemes/CentralDifferenceScheme.cpp
    src/Schemes/ConvectionSchemes/SecondOrderUpwindScheme.cpp
    src/Schemes/TimeSchemes/TimeScheme.cpp
    src/LinearSystem/Matrix.cpp
    src/LinearSystem/LinearSolvers.cpp
    src/LinearSystem/PETScRuntime.cpp
    src/Parallel/MPIRuntime.cpp
    src/Parallel/Comm.cpp
    src/Parallel/Reduce.cpp
    src/Parallel/GlobalIndex.cpp
    src/Parallel/MeshDecomposer.cpp
    src/Parallel/MeshDistributor.cpp
    src/Parallel/DecompositionChecker.cpp
    src/Solver/MomentumTransport.cpp
    src/Solver/Segregated.cpp
    src/Solver/SIMPLE.cpp
    src/Solver/PISO.cpp
    src/Models/Turbulence/TurbulenceModel.cpp
    src/Models/Turbulence/RANS.cpp
    src/Models/Turbulence/kOmegaSST.cpp
    src/PostProcessing/VTK/HDF5Wrapper.cpp
    src/PostProcessing/VTK/HDF5CellData.cpp
    src/PostProcessing/VTK/HDF5BoundaryData.cpp
    src/PostProcessing/PostProcess.cpp
    src/PostProcessing/DerivedFields.cpp
    src/PostProcessing/Forces.cpp
    src/Case/CaseReader.cpp
    src/Case/CaseConfiguration.cpp
    src/Case/RuntimeSelection.cpp
    src/Application/SolverSetup.cpp
    src/Application/CFDApplication.cpp
)

# --- Include directories (PUBLIC: propagated to the executable and tests) ---
target_include_directories(turblyze_core PUBLIC
    src/Mesh/
    src/Primitives/
    src/Fields/
    src/BoundaryConditions/
    src/BoundaryConditions/BoundaryTypes/
    src/Schemes/
    src/Schemes/Interpolation/
    src/Schemes/GradientSchemes/
    src/Schemes/ConvectionSchemes/
    src/Schemes/TimeSchemes/
    src/Models/
    src/Models/Turbulence/
    src/Case/
    src/Application/
    src/PostProcessing/
    src/PostProcessing/VTK/
    src/Solver/
    src/LinearSystem/
    src/Parallel/
)

# SYSTEM: /opt/homebrew/include also holds third-party headers whose
# warnings are not ours to fix
target_include_directories(turblyze_core SYSTEM PUBLIC
    ${METIS_INCLUDE_DIR}
)

turblyze_apply_build_flags(turblyze_core)

# PUBLIC so a test target linking turblyze_core builds with the matching
# Scalar/PetscScalar precision (a mismatch is a static_assert, not a silent bug)
if(TURBLYZE_USE_DOUBLE_PRECISION)
    target_compile_definitions(turblyze_core PUBLIC PROJECT_USE_DOUBLE_PRECISION)
endif()

# Link the solver core against HDF5, PETSc, MPI, METIS, and Eigen3 (PUBLIC so
# the executable and test targets inherit the same link line)
target_link_libraries(turblyze_core PUBLIC
    HDF5::HDF5
    PkgConfig::PETSC
    MPI::MPI_CXX
    ${METIS_LIBRARY}
    Eigen3::Eigen
)

# --- The solver executable: thin main() over the core library ---
add_executable(${PROJECT_NAME} src/main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE turblyze_core)

turblyze_apply_build_flags(${PROJECT_NAME})

# --- Unit and MPI tests (Catch2 + CTest). BUILD_TESTING defaults ON; configure
#     with -DBUILD_TESTING=OFF for a production-only build identical to before ---
include(CTest)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
