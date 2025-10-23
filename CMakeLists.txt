# Set the combiner name. The name must be the name of the subdirectory of the
# combiner.
set(COMBINER_NAME charm-fitter)

# Define all executables of this combiner. They have to be located in the main/
# subdirectory, and end with .cpp.
set(COMBINER_EXECUTABLES charmcombo ws-combo)

# -----------------------------------------------------------------------------
# Set directory names and includes
# -----------------------------------------------------------------------------

set(COMBINER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(COMBINER_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(COMBINER_MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/main)

include_directories(BEFORE ${ROOT_INCLUDE_DIR} ${CORE_INC_DIR}
                    ${COMBINER_INCLUDE_DIR})

# -----------------------------------------------------------------------------
# Build the library
# -----------------------------------------------------------------------------

set(COMBINER_LIB_SOURCES
    ${COMBINER_SOURCE_DIR}/CharmUtils.cpp
    ${COMBINER_SOURCE_DIR}/ParametersCharmCombo.cpp
    ${COMBINER_SOURCE_DIR}/PDF_AcpHH_LHCb_Run12.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_CLEO_K3pi_Kpipi0.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_Kpi_1d.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BES_Kpi.cpp
    ${COMBINER_SOURCE_DIR}/PDF_BinFlip.cpp
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

# -----------------------------------------------------------------------------
# Build the executables
# -----------------------------------------------------------------------------

set(COMBINER_LIBS ${COMBINER_LIB} ${CORE_LIB})
foreach(exec ${COMBINER_EXECUTABLES})
  add_executable(${exec} ${COMBINER_MAIN_DIR}/${exec}.cpp)
  target_link_libraries(${exec} ${COMBINER_LIBS})
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
