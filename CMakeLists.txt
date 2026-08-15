# Set the combiner name. The name must be the name of the subdirectory of the
# combiner.
set(COMBINER_NAME charm-fitter)

# Define all executables of this combiner. They have to be located in the main/
# subdirectory, and end with .cpp.
set(COMBINER_EXECUTABLES charm-combo ws-combo)

# -----------------------------------------------------------------------------
# Set directory names and includes
# -----------------------------------------------------------------------------

set(COMBINER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(COMBINER_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(COMBINER_MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/main)

# -----------------------------------------------------------------------------
# Build the library
# -----------------------------------------------------------------------------

set(COMBINER_LIB_SOURCES
    ${COMBINER_SOURCE_DIR}/CharmUtils.cpp
    ${COMBINER_SOURCE_DIR}/CharmParameters.cpp
    ${COMBINER_SOURCE_DIR}/PDF_AcpHH_LHCb_Run12.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_CLEO_K3pi_Kpipi0.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_Kpi_1d.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_Kpi_7d.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_Kpi.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BinFlip.cpp
    ${COMBINER_SOURCE_DIR}/PDF_Charm.cpp
    ${COMBINER_SOURCE_DIR}/PDF_Cleo.cpp
    ${COMBINER_SOURCE_DIR}/PDF_DY_pipipi0.cpp
    ${COMBINER_SOURCE_DIR}/PDF_DY_RS.cpp
    ${COMBINER_SOURCE_DIR}/PDF_DY.cpp
    ${COMBINER_SOURCE_DIR}/PDF_Fp_pipipi0.cpp
    ${COMBINER_SOURCE_DIR}/PDF_K3pi.cpp
    ${COMBINER_SOURCE_DIR}/PDF_Kpipi0.cpp
    ${COMBINER_SOURCE_DIR}/PDF_Kshh.cpp
    ${COMBINER_SOURCE_DIR}/PDF_RM.cpp
    ${COMBINER_SOURCE_DIR}/PDF_scan_DY_RS.cpp
    ${COMBINER_SOURCE_DIR}/PDF_WS_NoCPV.cpp
    ${COMBINER_SOURCE_DIR}/PDF_WS.cpp
    ${COMBINER_SOURCE_DIR}/PDF_XY.cpp
    ${COMBINER_SOURCE_DIR}/PDF_yCP_minus_yCP_KP.cpp
    ${COMBINER_SOURCE_DIR}/PDF_yCP_minus_yCP_RS.cpp
    ${COMBINER_SOURCE_DIR}/PDF_yCP_plus_yCP_RS.cpp
    ${COMBINER_SOURCE_DIR}/PDF_yCP.cpp)

set(COMBINER_LIB ${COMBINER_NAME}Components)
add_library(${COMBINER_LIB} SHARED ${COMBINER_LIB_SOURCES})
target_link_libraries(${COMBINER_LIB} ${CORE_LIB})
target_include_directories(${COMBINER_LIB} PUBLIC ${COMBINER_INCLUDE_DIR})

# -----------------------------------------------------------------------------
# Build the executables
# -----------------------------------------------------------------------------

foreach(exec ${COMBINER_EXECUTABLES})
  add_executable(${exec} ${COMBINER_MAIN_DIR}/${exec}.cpp)
  target_link_libraries(${exec} PUBLIC ${COMBINER_LIB})
endforeach()

# -----------------------------------------------------------------------------
# Install the binaries from the build directory back into the project
# subdirectory
# -----------------------------------------------------------------------------

foreach(exec ${COMBINER_EXECUTABLES})
  install(PROGRAMS $<TARGET_FILE:${exec}>
          DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endforeach()

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

foreach(folder ${OUTPUT_FOLDERS})
  install(DIRECTORY DESTINATION ${COMBINER_NAME}/${folder})
endforeach()

# -----------------------------------------------------------------------------
# Build the BLUE library and executables (optional)
# -----------------------------------------------------------------------------

option(
  BUILD_BLUE
  "Build the BLUE library and the executables for single-observable averages"
  ON)

if(BUILD_BLUE)
  set(BLUE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BLUE/include)
  set(BLUE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BLUE/src)
  set(BLUE_MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BLUE/main)

  set(BLUE_LIB Blue)
  add_library(${BLUE_LIB} SHARED ${BLUE_SOURCE_DIR}/Blue.cxx
                                 ${BLUE_SOURCE_DIR}/Utils.cpp)
  target_link_libraries(
    ${BLUE_LIB}
    PUBLIC ROOT::Core
           ROOT::Gpad
           ROOT::Graf
           ROOT::Hist
           ROOT::MathCore
           ROOT::Matrix
           ROOT::RIO)
  target_include_directories(${BLUE_LIB} PUBLIC ${BLUE_INCLUDE_DIR})

  # Silence compiler warnings for the BLUE library.
  if(MSVC)
    set_property(
      SOURCE ${BLUE_SOURCE_DIR}/Blue.cxx
      APPEND
      PROPERTY COMPILE_OPTIONS /w)
  else()
    set_property(
      SOURCE ${BLUE_SOURCE_DIR}/Blue.cxx
      APPEND
      PROPERTY COMPILE_OPTIONS -w)
  endif()

  root_generate_dictionary(G__${BLUE_LIB} BLUE/Blue.h MODULE ${BLUE_LIB}
                           LINKDEF ${BLUE_INCLUDE_DIR}/BLUE/LinkDef.h)

  set(BLUE_EXECUTABLES dy d0-to-ksks ycp ycp-minus-ycp-rs)
  foreach(exec ${BLUE_EXECUTABLES})
    add_executable(${exec} ${BLUE_MAIN_DIR}/${exec}.cpp)
    target_link_libraries(${exec} PUBLIC ${BLUE_LIB})
  endforeach()

  foreach(exec ${BLUE_EXECUTABLES})
    install(PROGRAMS $<TARGET_FILE:${exec}>
            DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin/BLUE)
  endforeach()
endif()
