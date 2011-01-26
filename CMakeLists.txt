cmake_minimum_required(VERSION 2.8.2)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

project(NanoITK)

include_directories(
  Code/BasicFilters
  Code/Common
)