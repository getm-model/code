add_library(output
            ascii_out.F90
#            calc_mean_fields.F90
#            diagnose.F90
#            diagnostic_variables.F90
            ncdf_out.F90
            nesting.F90
            output_processing.F90
            output.F90
           )
target_link_libraries(output PRIVATE domain 2d 3d field_manager output_manager futils)
