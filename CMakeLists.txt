#This file shows how to link to a MFEM (and PUMI) build
#installation using CMake
#it represents a simple 'CMakeLists.txt'
#file for a new project

cmake_minimum_required(VERSION 3.0.0)

project(inclusion_solver VERSION 1.0.0 LANGUAGES CXX)


# Let CMake know where to find custom FindFoo.cmake files
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

find_package(HYPRE REQUIRED)
find_package(METIS REQUIRED)


if (MFEM_PREFIX)
find_package(MFEM REQUIRED CONFIG PATHS ${MFEM_PREFIX} NO_DEFAULT_PATH)
else()
# IF MFEM_PREFIX was not specified, look in typical system directories,
# and also in CMAKE_PREFIX_PATH (environment variable)
find_package(
    MFEM   #package name, has to be MFEM
    REQUIRED #indicate that MFEM is really needed to compile
    CONFIG   #skip the 'MODULE' search system, save some time and confusion
    )
endif()

message(STATUS "MFEM INCLUDE is ${MFEM_INCLUDE_DIRS} ")
message(STATUS "MFEM LIBRARY is ${MFEM_LIBRARY_DIR} ")


add_library(inclusion_solver STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/inclusion_solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pfem_extras.cpp)
include_directories(${MFEM_INCLUDE_DIRS})
target_link_libraries(inclusion_solver ${MFEM_LIBRARY_DIR}/libmfem.a)

add_executable(inclusion ${CMAKE_CURRENT_SOURCE_DIR}/inclusion.cpp)
include_directories(${MFEM_INCLUDE_DIRS})
target_link_libraries(inclusion 
  ${MFEM_LIBRARY_DIR}/libmfem.a
  ${HYPRE_LIBRARIES}
  ${METIS_LIBRARIES}
  inclusion_solver)
