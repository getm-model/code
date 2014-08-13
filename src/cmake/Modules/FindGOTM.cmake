# Try to locate GOTM's installation prefix.
find_path(GOTM_PREFIX
  NAMES include/turbulence.mod
  PATHS "$ENV{LOCALAPPDATA}/gotm" "$ENV{APPDATA}/gotm" "$ENV{HOME}/local/gotm"
  DOC "Installation prefix for General Ocean Turbulence Models - gotm.net"
)

# Find GOTM turbulence library
find_library(GOTM_TURBULENCE NAMES turbulence
             HINTS ${GOTM_PREFIX}/lib
             DOC "GOTM turbulence library")

# Find GOTM utility library
find_library(GOTM_UTIL NAMES util
             HINTS ${GOTM_PREFIX}/lib
             DOC "GOTM uutility library")

set(GOTM_LIBRARIES "${GOTM_TURBULENCE} ${GOTM_UTIL}")

# Store configurable path of GOTM include directory
set(GOTM_INCLUDE_DIRS "${GOTM_PREFIX}/include"
    CACHE PATH
    "GOTM include directories")

mark_as_advanced(GOTM_LIBRARIES GOTM_INCLUDE_DIRS)

# Process default arguments (QUIET, REQUIRED)
include(FindPackageHandleStandardArgs) 
find_package_handle_standard_args (GOTM DEFAULT_MSG GOTM_LIBRARIES GOTM_INCLUDE_DIRS) 

# For backward compatibility:
#KBset(GOTM_LIBRARY GOTM_LIBRARIES)
#KBset(GOTM_INCLUDE_DIR GOTM_INCLUDE_DIRS)
