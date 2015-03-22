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
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")

# Set custom install folders here
set(PETSC_DIR "/home/manic/petsc_optimized")

# Double or float precision?
set(REAL "double")

# Use as ARRAY_ARGS ARRAY_SIZE while declaring arrays in functions arguments
# Must use -std=c99 if using the "restrict" keyword while compiling 
set(ARRAY_ARGS "const restrict static")

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

set(PROBLEM "linear_modes")
# Now choose problem specific options in the respective problem folder in
# src/problem/
# ------------------------------END OF USER OPTIONS---------------------------#

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Need to find all packages before starting compilation, so that all headers
# needed by the code are included.

# Should the external libraries needed by the internal libraries go into their
# own CMakeLists in their own folders?

# TODO: Only find packages depending on what options are set.

# System libraries
find_library(MATH_LIBRARIES m REQUIRED)

# External packages
find_package(PETSc  REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GSL    REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${PETSC_INCLUDES} ${GSL_INCLUDE_DIRS})

add_subdirectory(geometry)
add_subdirectory(physics)
add_subdirectory(reconstruct)
add_subdirectory(riemannsolver)
add_subdirectory(boundary)
add_subdirectory(gridzone)
add_subdirectory(timestepper)

set(PROBLEM_DIR ${CMAKE_SOURCE_DIR}/problem/${PROBLEM})
include(${PROBLEM_DIR}/CMakeLists.txt)
configure_file(${PROBLEM_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/inputs.h)
add_subdirectory(problem)

add_executable(grim grim.c grim.h)

target_link_libraries(grim geometry physics reconstruct boundary gridzone
                      timestepper
                      ${MATH_LIBRARIES} ${PETSC_LIBRARIES} ${GSL_LIBRARIES})

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("PETSc dir        : " ${PETSC_DIR})
message("")
