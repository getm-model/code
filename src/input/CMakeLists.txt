add_library(input OBJECT
            get_2d_bdy.F90
            get_2d_field.F90
            get_3d_bdy.F90
            get_3d_field.F90
            get_meteo_data.F90
            get_river_data.F90
            inquire_file.F90
            init_2d_bdy.F90
            init_3d_bdy.F90
            init_meteo_input.F90
            init_river_input.F90
            input.F90
            read_profile.F90
            read_topo_file.F90
           )
target_link_libraries(input PRIVATE ncdf 2d 3d field_manager netcdf)
