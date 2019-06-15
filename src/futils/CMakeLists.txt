configure_file("${CMAKE_CURRENT_SOURCE_DIR}/getm_version.F90.in" ${CMAKE_CURRENT_BINARY_DIR}/getm_version.F90)
if(GETM_EMBED_VERSION)
  add_custom_target(getm_version COMMAND
    ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_SOURCE_DIR}/getm_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/getm_version.F90 -P "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Retrieving description of last GETM commit..."
#KB    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/futils/getm_version.F90.in "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
    VERBATIM
  )
endif()

if(GETM_USE_PARALLEL)
   set(HALO_MPI halo_mpi.F90 ${MPI_Fortran_EXTRA_SOURCES})
endif(GETM_USE_PARALLEL)

add_library(futils OBJECT
#            git_revision.h
            c2x.F90
            check_3d_fields.F90
            cnv_2d.F90
            cnv_3d.F90
            col_interpol.F90
            eta_mask.F90
            exceptions.F90
            getm_error.F90
            getm_timers.F90
            grid_interpol.F90
            ${HALO_MPI}
            halo_zones.F90
            kbk_interpol.F90
            parallel.F90
            parameters.F90
            pos.F90
            read_par_setup.F90
            strip_string.F90
            time.F90
            to_2d_u.F90
            to_2d_vel.F90
            to_2d_v.F90
            to_3d_uu.F90
            to_3d_vel.F90
            to_3d_vv.F90
            tow.F90
            tridiagonal.F90
            ver_interpol.F90
            ${CMAKE_CURRENT_BINARY_DIR}/getm_version.F90
           )
if(GETM_EMBED_VERSION)
  add_dependencies(futils getm_version)
endif()
target_include_directories(futils PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules> $<INSTALL_INTERFACE:include>)
if(GETM_USE_PARALLEL)
  target_link_libraries(futils PRIVATE MPI::MPI_Fortran)
endif()
target_link_libraries(futils PRIVATE netcdf)
