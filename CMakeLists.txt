cmake_minimum_required(VERSION 2.8.2)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

project(NanoITK)

set(EXECUTABLE_OUTPUT_PATH "${NanoITK_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${NanoITK_BINARY_DIR}/lib")


#######################################
# Include ITK
#######################################

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
     "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


include_directories(
  Code/Algorithms
  Code/BasicFilters
  Code/Common
)

add_subdirectory( Tests )