set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# various compilation switches

if(GETM_USE_STATIC)
   add_definitions(-DSTATIC)
   configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/dimensions.h.in" ${CMAKE_CURRENT_BINARY_DIR}/dimensions.h)
   include_directories(${CMAKE_CURRENT_BINARY_DIR})
endif(GETM_USE_STATIC)

if(GETM_USE_DEBUG)
   add_definitions(-DDEBUG)
endif(GETM_USE_DEBUG)

if(GETM_USE_PROFILING)
   add_definitions(-DPROF)
endif(GETM_USE_PROFILING)

if(GETM_USE_FABM)
   add_definitions(-D_FABM_)
endif(GETM_USE_FABM)

if(GETM_USE_STRUCTURE_FRICTION)
   add_definitions(-DSTRUCTURE_FRICTION)
endif(GETM_USE_STRUCTURE_FRICTION)

if(GETM_FLAGS)
   add_definitions(${GETM_FLAGS})
endif(GETM_FLAGS)

add_definitions(-DFORTRAN95)

# used for naming executable
set(par_str "")
if(GETM_USE_PARALLEL)
   set(par_str "_parallel")
   add_definitions(-DGETM_PARALLEL)
endif(GETM_USE_PARALLEL)

# Specify coordinate type
if(NOT GETM_COORDINATE_TYPE)
   set (GETM_COORDINATE_TYPE "Cartesian" CACHE STRING "Choose coordinate type." FORCE)
   set_property(CACHE GETM_COORDINATE_TYPE PROPERTY STRINGS Cartesian Spherical Curvilinear)
endif()
if(${GETM_COORDINATE_TYPE} STREQUAL "Cartesian")
  add_definitions(-DCARTESIAN)
  set(exe_id "_cartesian${par_str}")
elseif(${GETM_COORDINATE_TYPE} STREQUAL "Spherical")
  add_definitions(-DSPHERICAL)
  set(exe_id "_spherical${par_str}")
elseif(${GETM_COORDINATE_TYPE} STREQUAL "Curvilinear")
  add_definitions(-DCURVILINEAR)
  set(exe_id "_curvilinear${par_str}")
endif()

# Customize compiler flags
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

# Global include directories
include_directories("${CMAKE_BINARY_DIR}"
                    "${PROJECT_SOURCE_DIR}/include"
                   )

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/2d)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/3d)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/domain)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/futils)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/getm)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/input)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/meteo)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/ncdf)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/output)

# Build GETM executable.
add_executable(getm_exe
               src/getm/main.F90
              )
set_property(TARGET getm_exe PROPERTY OUTPUT_NAME "getm${exe_id}")
target_link_libraries(getm_exe PRIVATE getm futils domain meteo 2d 3d ncdf input output output_manager gotm::turbulence)
if (NetCDF_STATIC_MSVC_BUILD)
  set_property(TARGET getm_exe PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmt\"")
endif()
if(GETM_USE_FABM)
   target_link_libraries(getm_exe PRIVATE gotm_fabm fabm)
endif()
if(GETM_USE_PARALLEL)
  target_link_libraries(getm_exe PRIVATE MPI::MPI_Fortran)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GETM in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GETM in")
endif()

# Install GETM executable
install(TARGETS getm_exe DESTINATION bin)
