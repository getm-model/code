cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Create main project.
project(getm Fortran)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# various compilation switches

option(USE_DEBUG "Debug compilation" OFF)
if(USE_DEBUG)
add_definitions(-DDEBUG)
endif(USE_DEBUG)

option(USE_PROFILING "Profiling compilation" OFF)
if(USE_PROFILING)
add_definitions(-DPROF)
endif(USE_PROFILING)

option(USE_STATIC "Static compilation" OFF)
if(USE_STATIC)
add_definitions(-DSTATIC)
endif(USE_STATIC)

option(USE_STRUCTURE_FRICTION "Include structure friction" OFF)
if(USE_STRUCTURE_FRICTION)
add_definitions(-DSTRUCTURE_FRICTION)
endif(USE_STRUCTURE_FRICTION)

add_definitions(-DFORTRAN95)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "-ffree-line-length-none")
  add_definitions(-DGFORTRAN)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  add_definitions(-DIFORT)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "-N 1023")
endif()

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
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
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GETM in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GETM in")
endif()

# Global include directories
include_directories("${PROJECT_SOURCE_DIR}/../include"
                    "${CMAKE_BINARY_DIR}/modules"
                   )
# Libraries added in alphabetic order - dependencies are dealt with below 
add_library(2d OBJECT
            2d/adv_arakawa_j7_2dh.F90
            2d/advection.F90
            2d/adv_fct_2dh.F90
            2d/adv_split_u.F90
            2d/adv_split_v.F90
            2d/adv_upstream_2dh.F90
            2d/bottom_friction.F90
            2d/cfl_check.F90
            2d/depth_update.F90
            2d/m2d.F90
            2d/momentum.F90
            2d/residual.F90
            2d/sealevel.F90
            2d/update_2d_bdy.F90
            2d/uv_advect.F90
            2d/uv_diff_2dh.F90
            2d/uv_diffusion.F90
            2d/variables_2d.F90
           )

add_library(3d OBJECT
            3d/adaptive_coordinates.F90
            3d/advection_3d.F90
            3d/adv_split_w.F90
            3d/bdy_3d.F90
            3d/bottom_friction_3d.F90
            3d/check_h.F90
            3d/coordinates.F90
            3d/eqstate.F90
            3d/general_coordinates.F90
#            3d/getm_bio.F90
#            3d/getm_fabm.F90
            3d/gotm.F90
            3d/hcc_check.F90
            3d/hybrid_coordinates.F90
            3d/internal_pressure.F90
            3d/ip_blumberg_mellor.F90
            3d/ip_blumberg_mellor_lin.F90
            3d/ip_chu_fan.F90
            3d/ip_shchepetkin_mcwilliams.F90
            3d/ip_song_wright.F90
            3d/ip_stelling_vankester.F90
            3d/ip_z_interpol.F90
            3d/m3d.F90
            3d/numerical_mixing.F90
            3d/physical_mixing.F90
            3d/preadapt_coordinates.F90
            3d/rivers.F90
            3d/salinity.F90
            3d/sigma_coordinates.F90
            3d/slow_bottom_friction.F90
            3d/slow_terms.F90
#            3d/spm.F90
            3d/ss_nn.F90
            3d/start_macro.F90
            3d/stop_macro.F90
            3d/stresses_3d.F90
#            3d/structure_friction_3d.F90
            3d/temperature.F90
            3d/tke_eps_advect_3d.F90
            3d/uu_momentum_3d.F90
            3d/uv_advect_3d.F90
            3d/uv_diffusion_3d.F90
            3d/variables_3d.F90
            3d/vv_momentum_3d.F90
            3d/ww_momentum_3d.F90
           )

add_library(domain OBJECT
            domain/bdy_spec.F90
            domain/domain.F90
            domain/have_bdy.F90
            domain/mirror_bdy_2d.F90
            domain/mirror_bdy_3d.F90
            domain/part_domain.F90
            domain/print_bdy.F90
            domain/uv_depths.F90
           )

add_library(futils OBJECT
            futils/c2x.F90
            futils/check_3d_fields.F90
            futils/cnv_2d.F90
            futils/cnv_3d.F90
            futils/col_interpol.F90
            futils/eta_mask.F90
            futils/exceptions.F90
            futils/getm_error.F90
            futils/getm_timers.F90
            futils/grid_interpol.F90
#            futils/halo_mpi.F90
            futils/halo_zones.F90
            futils/kbk_interpol.F90
            futils/parallel.F90
            futils/parameters.F90
            futils/pos.F90
            futils/read_par_setup.F90
            futils/time.F90
            futils/to_2d_u.F90
            futils/to_2d_vel.F90
            futils/to_2d_v.F90
            futils/to_3d_uu.F90
            futils/to_3d_vel.F90
            futils/to_3d_vv.F90
            futils/tow.F90
            futils/tridiagonal.F90
            futils/ver_interpol.F90
           )

add_library(getm OBJECT
            getm/cleanup.F90
            getm/initialise.F90
            getm/integration.F90
           )

add_library(input OBJECT
            input/get_2d_bdy.F90
            input/get_2d_field.F90
            input/get_3d_bdy.F90
            input/get_3d_field.F90
            input/get_meteo_data.F90
            input/get_river_data.F90
            input/init_2d_bdy.F90
            input/init_3d_bdy.F90
            input/init_meteo_input.F90
            input/init_river_input.F90
            input/input.F90
            input/read_profile.F90
            input/read_topo_file.F90
           )

add_library(meteo OBJECT
            meteo/exchange_coefficients.F90
            meteo/fluxes.F90
            meteo/meteo.F90
            meteo/short_wave_radiation.F90
           )

add_library(ncdf OBJECT
            ncdf/create_restart_ncdf.F90
            ncdf/grid_ncdf.F90
            ncdf/init_2d_ncdf.F90
            ncdf/init_3d_ncdf.F90
            ncdf/init_grid_ncdf.F90
            ncdf/init_mean_ncdf.F90
            ncdf/ncdf_2d_bdy.F90
            ncdf/ncdf_2d.F90
            ncdf/ncdf_3d_bdy.F90
            ncdf/ncdf_3d_bio_bdy.F90
            ncdf/ncdf_3d.F90
            ncdf/ncdf_close.F90
            ncdf/ncdf_common.F90
            ncdf/ncdf_get_field.F90
            ncdf/ncdf_in.F90
            ncdf/ncdf_mean.F90
            ncdf/ncdf_meteo.F90
            ncdf/ncdf_restart.F90
            ncdf/ncdf_rivers.F90
            ncdf/ncdf_topo.F90
            ncdf/open_restart_ncdf.F90
            ncdf/read_restart_ncdf.F90
            ncdf/save_2d_ncdf.F90
            ncdf/save_3d_ncdf.F90
            ncdf/save_grid_ncdf.F90
            ncdf/save_mean_ncdf.F90
            ncdf/set_attributes.F90
            ncdf/write_restart_ncdf.F90
           )

add_library(output OBJECT
            output/ascii_out.F90
            output/calc_mean_fields.F90
            output/diagnose.F90
            output/diagnostic_variables.F90
            output/ncdf_out.F90
            output/nesting.F90
            output/output.F90
           )

option(USE_ICE "Include ICE support" OFF)
if(USE_ICE)
  add_library(ice OBJECT
             )
endif(USE_ICE)

add_dependencies(domain futils)
add_dependencies(meteo futils domain)
add_dependencies(2d domain)
add_dependencies(3d meteo domain 2d)
add_dependencies(output futils domain 2d 3d meteo)
add_dependencies(ncdf futils domain output)
add_dependencies(input futils 2d 3d ncdf)
add_dependencies(getm input)

#KBset_property(TARGET util turbulence PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

find_package(NetCDF REQUIRED)
set_property(TARGET futils ncdf input APPEND PROPERTY INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIRS}")
add_definitions(-DNETCDF_FMT -DREAL_4B=real\(4\))

# Use GOTM_PREFIX to locate include and lib directories
  find_package(GOTM REQUIRED)
#  add_library(3d
#              3d/gotm.F90
#             )
  #KBadd_dependencies(gotm_fabm output)
#KB  target_link_libraries(gotm_fabm util)
    set_property(TARGET 3d getm APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTM_INCLUDE_DIRS}")
#KB  set_property(TARGET gotm_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTM_INCLUDE_DIRS}")
#KB  set_property(TARGET gotm_fabm PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

  #KBadd_dependencies(gotm gotm_fabm)
  #KBset_property(TARGET gotm APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")


#KBendif(USE_FABM)
# Use FABM_PREFIX to locate include and lib directories
#KBoption(USE_FABM "Include FABM support" ON)
#KBif(USE_FABM)
  #KBfind_package(FABM REQUIRED)
  #KBadd_library(gotm_fabm
              #KBfabm/gotm_fabm.F90
              #KBfabm/gotm_fabm_input.F90
              #KBfabm/gotm_fabm_output.F90
             #KB)
  #KBadd_dependencies(gotm_fabm output)
  #KBtarget_link_libraries(gotm_fabm util)
  #KBif(USE_NetCDF)
    #KBset_property(TARGET gotm_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIRS}")
  #KBendif(USE_NetCDF)
  #KBset_property(TARGET gotm_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
  #KBset_property(TARGET gotm_fabm PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
  #KBadd_definitions(-D_FABM_)

  #KBadd_dependencies(gotm gotm_fabm)
  #KBset_property(TARGET gotm APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
#KBendif(USE_FABM)

# Build GETM executable.
add_executable(getm_exe
               getm/main.F90
               $<TARGET_OBJECTS:getm>
               $<TARGET_OBJECTS:output>
               $<TARGET_OBJECTS:meteo>
               $<TARGET_OBJECTS:3d>
               $<TARGET_OBJECTS:2d>
               $<TARGET_OBJECTS:input>
               $<TARGET_OBJECTS:domain>
               $<TARGET_OBJECTS:ncdf>
               $<TARGET_OBJECTS:futils>
              )
set_property(TARGET getm_exe PROPERTY OUTPUT_NAME "getm")
target_link_libraries(getm_exe ${GOTM_TURBULENCE} ${GOTM_UTILITY} ${NetCDF_LIBRARIES})

#KBif(USE_NetCDF)
  #KBtarget_link_libraries(gotm_exe "${NetCDF_LIBRARIES}")
  #KBif (NetCDF_STATIC_MSVC_BUILD)
    #KBset_property(TARGET gotm_exe PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmtd\"")
    #KBtarget_compile_options(gotm_exe PRIVATE "/libs:static")
  #KBendif()
#KBendif(USE_NetCDF)

#KBif(USE_FABM)
  #KBtarget_link_libraries(gotm_exe gotm_fabm ${FABM_LIBRARIES})
#KBendif(USE_FABM)

# Install GETM library and include files, unless the host overruled this.
#KBinstall(TARGETS gotm_exe DESTINATION bin)
#KBinstall(TARGETS turbulence util DESTINATION lib)
#KBif(USE_FABM)
  #KBinstall(TARGETS gotm_fabm DESTINATION lib)
#KBendif()
#KBinstall(DIRECTORY ${CMAKE_BINARY_DIR}/modules/ DESTINATION include)

