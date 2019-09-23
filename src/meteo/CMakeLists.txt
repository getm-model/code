add_library(meteo
            exchange_coefficients.F90
            fluxes.F90
            meteo.F90
            solar_zenith_angle.F90
            short_wave_radiation.F90
            albedo_water.F90
           )
target_link_libraries(meteo PRIVATE domain futils)
