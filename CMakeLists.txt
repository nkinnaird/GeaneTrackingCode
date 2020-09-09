# ======================================================================
#  main build file
# ======================================================================

# ensure we are using a moden version of CMake
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

# project name - use all lowercase
PROJECT (gm2tracker)

# cetbuildtools contains our cmake modules
set( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
if( NOT CETBUILDTOOLS_VERSION )
  message( FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules" )
endif()

set( CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules ${CMAKE_MODULE_PATH} )
set( CMAKE_MODULE_PATH $ENV{CETBUILDTOOLS_DIR}/Modules ${CMAKE_MODULE_PATH} )

# include the necessary macros
include(CetCMakeEnv)
cet_cmake_env()


# g-2 wants the product name in the library name
set( ART_MAKE_PREPEND_PRODUCT_NAME TRUE)

# set compiler flags
cet_set_compiler_flags( DIAGS CAUTIOUS 
   WERROR
   EXTRA_FLAGS -pedantic
)

cet_report_compiler_flags()

# Locate some UPS areas (the versions are just minimum required) and setup include paths
# These are needed for Art stuff
find_ups_geant4(v4_10_3_p03agm2)
find_ups_product(cetbuildtools v7_00_03)
find_ups_product(art v2_10_03)
find_ups_product(fhiclcpp v4_06_05)
find_ups_product(messagefacility v2_01_06)
find_ups_product(gsl v2_4)
find_ups_product(eigen v3_3_4)
find_ups_product(artg4 v9_20_00)
find_ups_product(gm2dataproducts v9_20_00)
find_ups_product(gm2geom v9_20_00)
find_ups_product(gm2util v9_20_00)
find_ups_product(gm2trackerdaq v9_20_00)
find_ups_product(gm2unpackers v9_20_00)
find_ups_product(gm2midastoart v9_20_00)
find_ups_product(gm2aux v9_20_00)

cet_find_library( GSL NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )
cet_find_library( GSLCBLAS NAMES gslcblas PATHS ENV GSL_LIB NO_DEFAULT_PATH )

# macros for art_make and simple plugins
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# build these directories
add_subdirectory( analyses )
add_subdirectory( fitter )
add_subdirectory( filters )
add_subdirectory( utils )
add_subdirectory( reco )
add_subdirectory( readout )
add_subdirectory( cal )
add_subdirectory( teststand )
add_subdirectory( vizutils )
add_subdirectory( fcl )
add_subdirectory( ups )
add_subdirectory( test )
add_subdirectory( quality )

# packaging facility
include(UseCPack)
