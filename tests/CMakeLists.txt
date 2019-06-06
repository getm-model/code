# Begin testing section

include_directories(BEFORE ${PROJECT_SOURCE_DIR}/include)

# add_custom_target(test_all COMMENT "building tests - nothing yet :-)")
add_executable(test_speed_adv_upstream_2dh
               test_speed_adv_upstream_2dh.F90
              )

add_executable(test_rotation
               test_rotation.F90
               $<TARGET_OBJECTS:futils>
              )
target_link_libraries(test_rotation PRIVATE futils output_manager)
if (NetCDF_STATIC_MSVC_BUILD)
  set_property(TARGET test_rotation PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmt\"")
endif()
if(GETM_USE_PARALLEL)
  target_link_libraries(test_rotation PRIVATE MPI::MPI_Fortran)
endif()

add_custom_target(test_all WORKING_DIRECTORY tests COMMENT "building tests")
add_dependencies(test_all test_speed_adv_upstream_2dh test_rotation )

# End of testing section
