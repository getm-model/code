cmake_minimum_required(VERSION 3.0)

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
if(GOTM_PREFIX)
   message("Linking to pre-compiled GOTM libraries")
#   find_package(gotm REQUIRED HINTS "${GOTM_PREFIX}" "${DEFAULT_PREFIX}/../gotm")
   if(GETM_USE_FABM)
      find_package(fabm REQUIRED HINTS "${FABM_PREFIX}" "${GOTM_PREFIX}")
   endif()
   find_package(yaml REQUIRED HINTS "${GOTM_PREFIX}")
   find_package(flexout REQUIRED HINTS "${GOTM_PREFIX}")
   find_package(gotm REQUIRED HINTS "${GOTM_PREFIX}")
else()
   find_path(GOTM_BASE CMakeLists.txt HINTS ENV GOTM_BASE DOC "Path to GOTM source directory.")
   if (NOT GOTM_BASE)
      find_path(GOTM_BASE CMakeLists.txt PATHS ${CMAKE_CURRENT_LIST_DIR}/extern/gotm DOC "Path to GOTM source directory.")
      if (NOT GOTM_BASE)
         if (NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/extern/gotm/CMakeLists.txt")
            message(FATAL_ERROR "GOTM not found at extern/gotm. Please retrieve this submodule first by running \"git submodule update --init --recursive\" within your GETM source directory.")
         endif()
      endif()
   endif()
   message("Configuring GOTM from: " ${GOTM_BASE})
   add_subdirectory(${GOTM_BASE} extern/gotm EXCLUDE_FROM_ALL)
   set(GOTM_BUILD_LIBRARIES_ONLY ON)
   set(GOTM_USE_FABM ${GETM_USE_FABM} CACHE BOOL "Include support for Framework for Aquatic Biogeochemical Models (fabm.net)" FORCE)
   mark_as_advanced(GOTM_PREFIX GOTM_USE_FABM)
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

include(src/CMakeLists.txt)
add_subdirectory(tests EXCLUDE_FROM_ALL)
