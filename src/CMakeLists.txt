# To compile the code, make a seperate build folder and then do the following:
#
# cd /location/to/build
# cmake /location/to/grim/src
# make
#
# If a different compiler than the one cmake detects is needed, then the
# recommended way is to do the following:
#
# cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpic++ /location/to/grim/src 

cmake_minimum_required(VERSION 2.8)

project(grim)

# ------------------------------USER OPTIONS----------------------------------#
# Options for icc:
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3 -vec-threshold0 -qopt-report=5 -qopt-report-phase=vec -xhost")
# Options for gcc:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -g -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -g -fopenmp")
set(ARCH "OpenCL") # Choose CPU/OpenCL/CUDA

# Set custom install folders here
# 
# Configuration for SAVIO
#set(PETSC_DIR "/global/home/users/foucart/petsc_installation")
#set(YAML_ROOT_DIR "/global/home/users/manic/yaml-cpp/build")
#set(ArrayFire_ROOT_DIR "/global/scratch/manic/arrayfire/build_opencl_new/package")
#set(INTEL_MKL_ROOT_DIR "/global/software/sl-6.x86_64/modules/langs/intel/2016.1.150/mkl")

# Configuration for Stampede
#set(PETSC_DIR "/work/03669/tg829690/petsc_installation")
#set(YAML_ROOT_DIR "/global/home/users/manic/yaml-cpp/build")
#set(ArrayFire_ROOT_DIR "/work/03669/tg829690/arrayfire/build/package")

# Configuration for Deathstar
set(PETSC_DIR "/home/mani/Downloads/petsc_3.6.3_installation")
set(YAML_ROOT_DIR "/global/home/users/manic/yaml-cpp/build")
set(ArrayFire_ROOT_DIR "/home/mani/Downloads/arrayfire/build_gcc/package")
set(INTEL_MKL_ROOT_DIR "/opt/intel/mkl")

# Choose the problem to be solved:
#
# 1) linear_modes
# 2) shock_tests
# 3) advection_test
# 4) magnetized_field_loop_advection
# 5) orzag_tang
# 6) magnetized_explosion
# 7) atmosphere
# 8) bondi_inflow
# 10) anisotropic_conduction
# 11) torus
# 12) buoyancy_instabilities
# 13) bondi_viscous

set(PROBLEM "orzag_tang")
# Now choose problem specific options in the respective problem folder in
# src/problem/
# ------------------------------END OF USER OPTIONS---------------------------#

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(ArrayFire_DIR "${ArrayFire_ROOT_DIR}/share/ArrayFire/cmake")

# Needed to link python interfaces to precompiled .so libraries
set(BUILD_SHARED_LIBS "ON")

# Need to find all packages before starting compilation, so that all headers
# needed by the code are included.

# Should the external libraries needed by the internal libraries go into their
# own CMakeLists in their own folders?

# TODO: Only find packages depending on what options are set.

# System libraries
find_library(MATH_LIBRARIES m REQUIRED)

# External packages
find_package(PETSc REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)
find_package(YAML REQUIRED)
find_package(ArrayFire REQUIRED)
find_package(LAPACKE REQUIRED)

if (ARCH STREQUAL "CPU")
  set(ArrayFire_LIBRARIES ${ArrayFire_CPU_LIBRARIES})
elseif (ARCH STREQUAL "CUDA")
  set(ArrayFire_LIBRARIES ${ArrayFire_CUDA_LIBRARIES})
elseif (ARCH STREQUAL "OpenCL")
  set(ArrayFire_LIBRARIES ${ArrayFire_OpenCL_LIBRARIES})
endif()

include(UseCython)
include_directories(${PETSC_INCLUDES} 
                    ${YAML_INCLUDES}
                    ${ArrayFire_INCLUDE_DIRS}
                    ${PYTHON_NUMPY_INCLUDE_DIR}
                    ${LAPACKE_INCLUDES}
                    grid
                    geometry
                    reconstruction
                    boundary
                    physics
                    timestepper
                   )

add_subdirectory(grid)
add_subdirectory(geometry)
add_subdirectory(physics)
add_subdirectory(reconstruction)
add_subdirectory(boundary)
add_subdirectory(timestepper)

set(PROBLEM_DIR ${CMAKE_SOURCE_DIR}/problem/${PROBLEM})
add_subdirectory(problem)

add_executable(grim grim.cpp grim.hpp params.hpp)

target_link_libraries(grim grid geometry physics reconstruction
                      timestepper problem boundary
                      ${MATH_LIBRARIES} 
                      ${PETSC_LIBRARIES}
                      ${YAML_LIBRARIES}
                      ${ArrayFire_LIBRARIES}
                      ${LAPACK_LIBRARIES}
                     )

set(NUM_PROCS 4)
set(N1_test   32)
set(N2_test   32)
set(N3_test   32)
enable_testing()
# 1D
add_test(communication_1D_X1_left_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k X1_left
        )
add_test(communication_1D_X1_right_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k X1_right
        )
# 2D
add_test(communication_2D_X1_left_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k X1_left
        )
add_test(communication_2D_X1_right_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k X1_right
        )
add_test(communication_2D_X2_top_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k X2_top
        )
add_test(communication_2D_X2_bottom_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k X2_bottom
        )
# 3D
add_test(communication_3D_X1_left_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X1_left
        )
add_test(communication_3D_X1_right_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X1_right
        )
add_test(communication_3D_X2_top_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X2_top
        )
add_test(communication_3D_X2_bottom_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X2_bottom
        )
add_test(communication_3D_X3_front_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X3_front
        )
add_test(communication_3D_X3_back_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X3_back
        )

add_test(X1Coords_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k X1Coords
        )

add_test(X1Coords_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k X1Coords
        )

add_test(X2Coords_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k X2Coords
        )

add_test(X1Coords_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X1Coords
        )

add_test(X2Coords_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X2Coords
        )

add_test(X3Coords_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/grid/test_mpi.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k X3Coords
        )

# geometry tests
add_test(Minkowski_params_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_params
        )
add_test(Minkowski_gCov_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_gCov
        )
add_test(Minkowski_gCon_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_gCon
        )
add_test(Minkowski_g_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_g
        )
add_test(Minkowski_alpha_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_alpha
        )
add_test(Minkowski_params_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_params
        )
add_test(Minkowski_gCov_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_gCov
        )
add_test(Minkowski_gCon_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_gCon
        )
add_test(Minkowski_g_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_g
        )
add_test(Minkowski_alpha_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_alpha
        )
add_test(Minkowski_params_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_params
        )
add_test(Minkowski_gCov_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_gCov
        )
add_test(Minkowski_gCon_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_gCon
        )
add_test(Minkowski_g_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_g
        )
add_test(Minkowski_alpha_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k minkowski_alpha
        )
# geometry: Kerr metric
set(H_SLOPE         0.3)
set(BLACK_HOLE_SPIN 0.9375)
add_test(Kerr_params_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_params
        )
add_test(Kerr_gCov_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gCov
        )
add_test(Kerr_gCon_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gCon
        )
add_test(Kerr_g_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_g
        )
add_test(Kerr_alpha_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_alpha
        )
add_test(Kerr_connection_coeffs_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gammaUpDownDown
        )
add_test(Kerr_params_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_params
        )
add_test(Kerr_gCov_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gCov
        )
add_test(Kerr_gCon_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gCon
        )
add_test(Kerr_g_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_g
        )
add_test(Kerr_alpha_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_alpha
        )
add_test(Kerr_connection_coeffs_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gammaUpDownDown
        )
add_test(Kerr_params_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_params
        )
add_test(Kerr_gCov_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gCov
        )
add_test(Kerr_gCon_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gCon
        )
add_test(Kerr_g_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_g
        )
add_test(Kerr_alpha_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_alpha
        )
add_test(Kerr_connection_coeffs_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/geometry/test_geometry.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --hSlope=${H_SLOPE} --blackHoleSpin=${BLACK_HOLE_SPIN}
         --build_path=${CMAKE_BINARY_DIR} -k modifiedKerrSchild_gammaUpDownDown
        )

# boundary tests
add_test(X1_left_outflow_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X1Left
        )
add_test(X1_left_mirror_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X1Left
        )
add_test(X1_right_outflow_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X1Right
        )
add_test(X1_right_mirror_1D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=1
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X1Right
        )
# 2D
add_test(X1_left_outflow_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X1Left
        )
add_test(X1_left_mirror_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X1Left
        )
add_test(X1_right_outflow_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X1Right
        )
add_test(X1_right_mirror_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X1Right
        )
add_test(X2_bottom_outflow_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X2Bottom
        )
add_test(X2_bottom_mirror_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X2Bottom
        )
add_test(X2_top_outflow_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X2Top
        )
add_test(X2_top_mirror_2D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=2
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X2Top
        )
# 3D
add_test(X1_left_outflow_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X1Left
        )
add_test(X1_left_mirror_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X1Left
        )
add_test(X1_right_outflow_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X1Right
        )
add_test(X1_right_mirror_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X1Right
        )
add_test(X2_bottom_outflow_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X2Bottom
        )
add_test(X2_bottom_mirror_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X2Bottom
        )
add_test(X2_top_outflow_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X2Top
        )
add_test(X2_top_mirror_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X2Top
        )
add_test(X3_back_outflow_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X3Back
        )
add_test(X3_back_mirror_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X3Back
        )
add_test(X3_front_outflow_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k outflow_X3Front
        )
add_test(X3_front_mirror_3D_${NUM_PROCS}_procs
         mpirun -np ${NUM_PROCS} 
         py.test  ${CMAKE_SOURCE_DIR}/boundary/test_boundary.py
         --N1=${N1_test} --N2=${N2_test} --N3=${N3_test} --dim=3
         --build_path=${CMAKE_BINARY_DIR} -k mirror_X3Front
        )

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("Problem          : " ${PROBLEM})
message("Architecture     : " ${ARCH})
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("CXX Compiler     : " ${CMAKE_CXX_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("CXX_FLAGS        : " ${CMAKE_CXX_FLAGS})
message("PETSc dir        : " ${PETSC_DIR})
message("ArrayFire dir    : " ${ArrayFire_ROOT_DIR})
message("MKL dir          : " ${LAPACKE_INCLUDES})
message("NumPy dir        : " ${PYTHON_NUMPY_INCLUDE_DIR})
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("")
