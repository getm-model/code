@rem Script to configure the getm executable uing CMake

@set old=%cd%
@rem echo %old%

@echo Build directory:
@if "%build_dir%"=="" ( @set build_dir=%TEMP%\build\getm ) else ( @echo build_dir is set )
@echo %build_dir%
@if not EXIST "%build_dir%\." ( @mkdir "%build_dir%" )
@chdir "%build_dir%"

@echo Base directories:
@set GETM_BASE=%USERPROFILE%\source\repos\\GETM\code
@echo %GETM_BASE%

@echo Default Fortran compiler is ifort
@set compiler=ifort

@echo Install directory:
@set install_prefix=%APPDATA%\getm
@echo %install_prefix%

@echo Ready to configure:
FOR %%c IN (Cartesian Spherical Curvilinear) DO (
   @IF NOT EXIST "%%c\." ( mkdir "%%c" )
   @chdir "%%c"
   cmake "%GETM_BASE%" ^
         -DGETM_EMBED_VERSION=on ^
         -DGETM_USE_FABM=on ^
         -DCMAKE_Fortran_COMPILER=%compiler% ^
         -DGETM_USE_PARALLEL=off ^
         -DGETM_COORDINATE_TYPE=%%c ^
         -DCMAKE_INSTALL_PREFIX="%install_prefix%"

   @chdir ..\
)

@pause

@chdir %old%
