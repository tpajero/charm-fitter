cmake_minimum_required(VERSION 3.14)

# The name must be the name of the subdirectory of the combiner.
set(COMBINER_NAME charmcombo)

# Define all executables of this combiner.
#  They have to be located in the main/ subdirectory, and
#  end with .cpp.
set(COMBINER_EXECUTABLES
 	${COMBINER_NAME}_phenomenological
	${COMBINER_NAME}_theoretical
	${COMBINER_NAME}_superweak
)

######################################
#
# set directory names and includes
#
######################################

set(COMBINER_INCLUDE_DIR ./include)
set(COMBINER_SOURCE_DIR ./src)
set(COMBINER_MAIN_DIR ./main)

include_directories(
	BEFORE 
	${ROOT_INCLUDE_DIRS}
	${CORE_INC_DIR}
	${COMBINER_INCLUDE_DIR}
)

######################################
#
# build the library
#
######################################

aux_source_directory(${COMBINER_SOURCE_DIR} COMBINER_LIB_SOURCES)
add_library(${COMBINER_NAME}Components SHARED ${COMBINER_LIB_SOURCES})
target_link_libraries(${COMBINER_NAME}Components ${ROOT_LIBRARIES} ${PROJECT_NAME}Components)

######################################
#
# build the executables
#
######################################

set(COMBINER_LIBS
	${COMBINER_NAME}Components  
	${PROJECT_NAME}Components
)
foreach(exec ${COMBINER_EXECUTABLES})
	add_executable(${exec} ${COMBINER_MAIN_DIR}/${exec}.cpp)
	target_link_libraries(${exec} ${COMBINER_LIBS})
endforeach()

######################################
#
# install the binaries from the build directory back into the project subdirectory
#
######################################

foreach(exec ${COMBINER_EXECUTABLES})
	install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/${exec}
          DESTINATION ${COMBINER_NAME}/bin/)
endforeach()

######################################
#
# create output directories
#
######################################

install(DIRECTORY DESTINATION ${COMBINER_NAME}/root         )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots        )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/C      )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/cl     )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/dot    )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/eps    )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/par    )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/pdf    )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/png    )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/root   )
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/scanner)

