#   Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(Optimization)

# Load UserSettings.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
 set(BUILD_STYLE "standalone")
 include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
 set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
 include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
 include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
 STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
 set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
 # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
 set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
 set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIBROOT)
 set(LIBROOT "${CODEROOT}/lib")
endif()
if(NOT BINROOT)
 set(BINROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
 set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
 install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
 set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
 set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
 set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
 get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
 add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Include the top-level directories.
include_directories(AFTER
 "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
 include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
 set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
 set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
 set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.55.0
            COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
 include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
 include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

 # Find CSPICE library on local system.
 find_package(Spice)

 # Include CSpice directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
 endif( )

if(NOT USE_NRLMSISE00)
 message(STATUS "NRLMSISE-00 disabled!")
 add_definitions(-DUSE_NRLMSISE00=0)
else()
 message(STATUS "NRLMSISE-00 enabled!")
 add_definitions(-DUSE_NRLMSISE00=1)
 # Find USE_NRLMSISE00 library on local system.
 find_package(NRLMSISE00)

 # Include NRLMSISE00 directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
 endif( )
endif( )

# Find PaGMO library on local system.
if( USE_PAGMO )
 find_package(PaGMO2)

 # Include PaGMO directories.
 if(NOT APPLE OR APPLE_INCLUDE_FORCE)
   include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
 else( )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
 endif( )
endif( )

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(tudatLinkLibraries)

## Add problem definition folder
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Problems")


#add_subdirectory("${PROJECT_SOURCE_DIR}/customlibs/fftw")
#add_library(fftw3 STATIC "${PROJECT_SOURCE_DIR}/customlibs/fftw/fftw3.h")

# Add default model settings
add_executable(T_defaultOpt "${SRCROOT}/T_DefaultOpt.cpp")
setup_executable_target(T_defaultOpt "${SRCROOT}")
target_link_libraries(T_defaultOpt pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Add default model settings
add_executable(T_costFunctionTransferTester "${SRCROOT}/T_CostFunctionTransferTester.cpp")
setup_executable_target(T_costFunctionTransferTester "${SRCROOT}")
target_link_libraries(T_costFunctionTransferTester fftw3 swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Comparison of optimization algorithms
add_executable(T_OptComparison "${SRCROOT}/T_OptComparison.cpp")
setup_executable_target(T_OptComparison "${SRCROOT}")
target_link_libraries(T_OptComparison pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Comparison of optimization algorithms
add_executable(T_FineGrainMethods "${SRCROOT}/T_OptFineGrainMethods.cpp")
setup_executable_target(T_FineGrainMethods "${SRCROOT}")
target_link_libraries(T_FineGrainMethods pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Comparison of optimization algorithms
add_executable(T_CoarseGrainMethods "${SRCROOT}/T_OptCoarseGrainMethods.cpp")
setup_executable_target(T_CoarseGrainMethods "${SRCROOT}")
target_link_libraries(T_CoarseGrainMethods pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Comparison of using multiple different optimization algorithms
add_executable(T_CoarseGrainMultiMethods "${SRCROOT}/T_OptCoarseGrainMultiMethods.cpp")
setup_executable_target(T_CoarseGrainMultiMethods "${SRCROOT}")
target_link_libraries(T_CoarseGrainMultiMethods pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# General algorithm to run on the Eudoxos server
add_executable(T_EudoxosOptimisation "${SRCROOT}/T_EudoxosOptimisation.cpp")
setup_executable_target(T_EudoxosOptimisation "${SRCROOT}")
target_link_libraries(T_EudoxosOptimisation pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Propagate the results of a optimisation datafile for a certain duration
add_executable(T_PropagateOptimisation "${SRCROOT}/T_PropagateOptimisationResult.cpp")
setup_executable_target(T_PropagateOptimisation "${SRCROOT}")
target_link_libraries(T_PropagateOptimisation ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})

# Second order optimisation, optimise relative velocities for long-term stability.
add_executable(T_SecondOrderOptimisation "${SRCROOT}/T_SecondOrderOptimization.cpp")
setup_executable_target(T_SecondOrderOptimisation "${SRCROOT}")
target_link_libraries(T_SecondOrderOptimisation pagmo swarm_optimization_library ${TUDAT_PROPAGATION_LIBRARIES}
    pthread ${Boost_LIBRARIES})
