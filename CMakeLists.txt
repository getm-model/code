cmake_minimum_required(VERSION 3.10)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Create main project.
project(getm VERSION 2.5.0 LANGUAGES Fortran)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Set default installation prefix.
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/getm")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/getm")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/getm")
endif()

option(GETM_USE_FABM "Include support for Framework for Aquatic Biogeochemical Models (fabm.net)" OFF)
set(GOTM_EMBED_VERSION ${GETM_EMBED_VERSION} CACHE BOOL "Embed GOTM version information" FORCE)
find_path(GOTM_BASE src/gotm/gotm.F90 DOC "Path to GOTM source directory.")
if(GOTM_BASE)
  message("Building GOTM related libraries")
  set(GOTM_BUILD_LIBRARIES_ONLY ON)
  set(GOTM_USE_FABM ${GETM_USE_FABM} CACHE BOOL "Include support for Framework for Aquatic Biogeochemical Models (fabm.net)" FORCE)
  add_subdirectory(${GOTM_BASE} ${CMAKE_BINARY_DIR}/extern/gotm)
  mark_as_advanced(GOTM_PREFIX GOTM_USE_FABM)
else()
  message("Linking to pre-compiled GOTM libraries")
  find_package(yaml REQUIRED HINTS "${GOTM_PREFIX}" "${DEFAULT_PREFIX}/../gotm")
  find_package(flexout REQUIRED HINTS "${GOTM_PREFIX}" "${DEFAULT_PREFIX}/../gotm")
  find_package(gotm REQUIRED HINTS "${GOTM_PREFIX}" "${DEFAULT_PREFIX}/../gotm")
endif()

option(GETM_USE_PARALLEL "Enable parallel computation" OFF)
if(GETM_USE_PARALLEL)
   find_package(MPI REQUIRED)
endif(GETM_USE_PARALLEL)

# Various compilation options
option(GETM_USE_STATIC "Static compilation (domain size is set at compile time)" OFF)

option(GETM_USE_DEBUG "Enable debug output" OFF)

option(GETM_USE_PROFILING "Profiling compilation" OFF)

option(GETM_USE_STRUCTURE_FRICTION "Include structure friction" OFF)

# Customize compiler flags
#[[---
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  add_definitions(-DGFORTRAN)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  add_definitions(-DIFORT)
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
    set(CMAKE_Fortran_STACK_SIZE 128000000)
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
endif()
---]]

add_subdirectory(src)
add_subdirectory(tests EXCLUDE_FROM_ALL)
