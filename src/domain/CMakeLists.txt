add_library(domain
            bdy_spec.F90
            domain.F90
            have_bdy.F90
            mirror_bdy_2d.F90
            mirror_bdy_3d.F90
            part_domain.F90
            print_bdy.F90
            uv_depths.F90
           )
target_link_libraries(domain PRIVATE futils input)
