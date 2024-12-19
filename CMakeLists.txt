cmake_minimum_required(VERSION 3.14)

project(CharmFitter VERSION 1.0.0)

# Install destinations can be passed via the command line:
#  cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=.... -DINCLUDE_OUTPUT_DIRECTORY:PATH=....
#  Otherwise, set default install destinations
if(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
  set(LIBRARY_OUTPUT_DIRECTORY ./lib)
endif()
if(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
  set(INCLUDE_OUTPUT_DIRECTORY ./include)
endif()

# Include cmake macros and settings
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(settings)
include(macros)
include(combiners)

# Use ccache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# Find FreeType
find_package(Freetype REQUIRED)

# Find ROOT
#  ROOT v6.20 ensures that the compiler is at least C++11
set(ROOT_ADDITIONAL_LIBRARIES RooFit RooFitCore RooFitMore Html Minuit Thread RooStats Gui TreePlayer GenVector)
find_package(ROOT 6.20 CONFIG REQUIRED COMPONENTS ${ROOT_ADDITIONAL_LIBRARIES})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O3 ${ROOT_CXX_FLAGS}")

# Find Boost
set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "Enable fix for FindBoost.cmake")
find_package(Boost REQUIRED)
if(Boost_FOUND)
  message(STATUS "Boost is found at: ${Boost_INCLUDE_DIR}")
  include_directories(${Boost_INCLUDE_DIR})
endif()

# Find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  # Add combiner modules to Doxygen
  set(COMBINER_MODULES_DOXYGEN_SOURCES)
  foreach(dir ${COMBINER_MODULES})
    set(COMBINER_MODULES_DOXYGEN_SOURCES "${COMBINER_MODULES_DOXYGEN_SOURCES}
      INPUT += ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include
      INPUT += ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/src"
    )
  endforeach()
  # Configure Doxygen config file
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  # Add target
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
else()
  message(STATUS "Doxygen not installed or found -- API documentation can't be generated!")
endif(DOXYGEN_FOUND)

######################################
#
# set includes
#
######################################

set(CORE_INC_DIR ./core/include)

include_directories(
  BEFORE
  ${CORE_INC_DIR}
  ${ROOT_INCLUDE_DIRS}
)

link_directories(${ROOT_LIBRARY_DIR})

aux_source_directory(./core/src CORE_LIB_SOURCES)
list(REMOVE_ITEM CORE_LIB_SOURCES "./core/src/SharedArray.cpp")

######################################
#
# generate ROOT dictionaries
#
######################################

set(CORE_DICTIONARY_SOURCES
  # RooBinned2DBicubicBase.h
  RooCrossCorPdf.h
  RooHistPdfAngleVar.h
  RooHistPdfVar.h
  RooPoly3Var.h
  RooPoly4Var.h
  RooSlimFitResult.h
  RooMultiPdf.h
)
set(CORE_DICTIONARY_CPP_SOURCES)
foreach(file ${CORE_DICTIONARY_SOURCES})
  cmake_path(REPLACE_EXTENSION file ".cpp")
  if ((NOT ${file} EQUAL "RooBinned2DBicubicBase.cpp") AND (NOT ${file} EQUAL "SharedArray.cpp"))
    set(CORE_DICTIONARY_CPP_SOURCES "./core/src/${file}" ${CORE_DICTIONARY_CPP_SOURCES})
  endif()
endforeach()
message("CPP sources: ${CORE_DICTIONARY_CPP_SOURCES}")

add_library(${PROJECT_NAME}Dict ${CORE_DICTIONARY_CPP_SOURCES} G__${PROJECT_NAME}Dict.cxx)

# Attach dictionaries to the executable. First, tell it where to look for headers required by the dictionaries:
target_include_directories(${PROJECT_NAME}Dict PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
# Then generate dictionaries and add them as a dependency of the executable (via the MODULE parameter):
root_generate_dictionary(G__${PROJECT_NAME}Dict
  ${CORE_DICTIONARY_SOURCES}
  MODULE ${PROJECT_NAME}Dict
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/core/include/coreLinkDef.h
)

######################################
#
# build core library
#
######################################

# Add the dictionary to the list of library sources
set(CORE_LIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/G__${PROJECT_NAME}Dict.cxx ${CORE_LIB_SOURCES})

add_shared_library(${PROJECT_NAME}Components ${CORE_LIB_SOURCES})
add_dependencies(${PROJECT_NAME}Components ${PROJECT_NAME}Dict)
target_link_libraries(${PROJECT_NAME}Components ${ROOT_LIBRARIES} ${FREETYPE_LIBRARIES})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}Dict_rdict.pcm
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lib)

######################################
#
# subdirectorys for the combiners
#
######################################

foreach(dir ${COMBINER_MODULES})
  add_subdirectory(${dir})
endforeach()
