# Try to locate GOTM's installation prefix.
find_path(GOTM_PREFIX
  NAMES include/turbulence.mod
  HINTS "$ENV{GOTM_PREFIX}"
  PATHS "$ENV{LOCALAPPDATA}/gotm" "$ENV{APPDATA}/gotm" "$ENV{HOME}/local/gotm"
  DOC "Installation prefix for General Ocean Turbulence Models - gotm.net"
)

# Find GOTM/FABM coupling library if USE_FABM
if(GETM_USE_FABM)
  find_library(GOTM_FABM NAMES gotm_fabm
               HINTS ${GOTM_PREFIX}/lib
               DOC "GOTM-FABM library")
endif(GETM_USE_FABM)

# Find GOTM turbulence library
find_library(GOTM_TURBULENCE NAMES turbulence
             HINTS ${GOTM_PREFIX}/lib
             DOC "GOTM turbulence library")

# Find GOTM utility library
find_library(GOTM_UTIL NAMES util
             HINTS ${GOTM_PREFIX}/lib
             DOC "GOTM utility library")

set(GOTM_LIBRARIES ${GOTM_FABM} ${GOTM_TURBULENCE} ${GOTM_UTIL})

# Store configurable path of GOTM include directory
find_path(GOTM_INCLUDE_DIRS
          NAMES turbulence.mod
          HINTS ${GOTM_PREFIX}/include
          DOC "GOTM include directories"
)

mark_as_advanced(GOTM_LIBRARIES GOTM_INCLUDE_DIRS GOTM_TURBULENCE GOTM_UTIL GOTM_FABM)

# Process default arguments (QUIET, REQUIRED)
include(FindPackageHandleStandardArgs) 
find_package_handle_standard_args (GOTM DEFAULT_MSG GOTM_LIBRARIES GOTM_INCLUDE_DIRS) 
