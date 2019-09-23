add_library(2d
            adv_arakawa_j7_2dh.F90
            advection.F90
            adv_fct_2dh.F90
            adv_split_u.F90
            adv_split_v.F90
            adv_upstream_2dh.F90
            bottom_friction.F90
            cfl_check.F90
            depth_update.F90
            m2d.F90
            momentum.F90
            residual.F90
            sealevel.F90
            update_2d_bdy.F90
            uv_advect.F90
            uv_diff_2dh.F90
            uv_diffusion.F90
            variables_2d.F90
           )
target_compile_definitions(2d PRIVATE REAL_4B=real\(4\))
target_link_libraries(2d PRIVATE domain field_manager)
