add_library(getm OBJECT
            cleanup.F90
            compilation_options.F90
            initialise.F90
            integration.F90
            register_all_variables.F90
            print_version.F90
           )
target_link_libraries(getm PRIVATE futils domain input output_manager gotm::util)
if (GETM_USE_FABM)
   target_link_libraries(getm PRIVATE fabm)
endif()
