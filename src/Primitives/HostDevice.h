/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HostDevice.h
 * @brief HD execution-space macro for functions shared by CPU and GPU code
 *
 * @details Functions marked HD compile for both the CPU (host) and the GPU
 * (device) when the translation unit is processed by nvcc, which defines
 * __CUDACC__ and requires the __host__ __device__ qualifiers. Under any
 * other compiler HD expands to nothing, so CPU-only builds are unchanged.
 *
 * Function bodies that must differ between the two targets (e.g. error
 * handling — device code cannot throw or reach ErrorHandler) branch on
 * __CUDA_ARCH__, which is defined only during nvcc's device compilation
 * pass; the host path is preserved verbatim.
 *****************************************************************************/

#pragma once

// ************************* Execution-Space Macro ****************************

#ifdef __CUDACC__
    #define HD __host__ __device__
#else
    #define HD
#endif
