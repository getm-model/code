if(GETM_USE_FABM)
   set (GETM_FABM getm_fabm.F90)
   set(FABM_EMBED_VERSION ${GETM_EMBED_VERSION} CACHE BOOL "Embed FABM version information" FORCE)
endif(GETM_USE_FABM)

if(GETM_USE_STRUCTURE_FRICTION)
   set (STRUCTURE_FRICTION structure_friction_3d.F90)
endif(GETM_USE_STRUCTURE_FRICTION)

add_library(3d OBJECT
            adaptive_coordinates.F90
            advection_3d.F90
            adv_split_w.F90
            bdy_3d.F90
            bottom_friction_3d.F90
            check_h.F90
            coordinates.F90
            eqstate.F90
            general_coordinates.F90
#            getm_bio.F90
            ${GETM_FABM}
            gotm.F90
            hcc_check.F90
            hybrid_coordinates.F90
            internal_pressure.F90
            ip_blumberg_mellor.F90
            ip_blumberg_mellor_lin.F90
            ip_chu_fan.F90
            ip_shchepetkin_mcwilliams.F90
            ip_song_wright.F90
            ip_stelling_vankester.F90
            ip_z_interpol.F90
            m3d.F90
            numerical_mixing.F90
            physical_mixing.F90
            preadapt_coordinates.F90
            rivers.F90
            salinity.F90
            sigma_coordinates.F90
            slow_bottom_friction.F90
            slow_terms.F90
#            spm.F90
            ss_nn.F90
            start_macro.F90
            stop_macro.F90
            stresses_3d.F90
            ${STRUCTURE_FRICTION}
            temperature.F90
            tke_eps_advect_3d.F90
            uu_momentum_3d.F90
            uv_advect_3d.F90
            uv_diffusion_3d.F90
            variables_3d.F90
            vv_momentum_3d.F90
            ww_momentum_3d.F90
           )
target_link_libraries(3d PRIVATE domain meteo 2d field_manager futils gotm::turbulence)
if (GETM_USE_FABM)
   target_link_libraries(3d PRIVATE fabm)
endif()
