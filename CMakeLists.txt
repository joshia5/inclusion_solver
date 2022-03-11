#This file shows how to link to a MFEM (and PUMI) build
#installation using CMake
#it represents a simple 'CMakeLists.txt'
#file for a new project

cmake_minimum_required(VERSION 3.0.0)

project(inclusion_solver VERSION 1.0.0 LANGUAGES CXX)


# Let CMake know where to find custom FindFoo.cmake files
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

if(USE_OMEGAH)
  find_package(Omega_h REQUIRED)
endif()
if(NOT MFEM_USE_CUDA)
  #HACK: Clear the omegah compilation flags that it passes to cuda.
  if(USE_OMEGAH)
    set_property(TARGET Omega_h::omega_h PROPERTY INTERFACE_COMPILE_OPTIONS "")
  endif()
endif()

if(MFEM_USE_CUDA)
  if(NOT OMEGA_H_USE_CUDA)
    #TODO we may want to support mfem+cuda and omegah without cuda at some
    # point... possibly for debugging
    message(FATAL_ERROR "MFEM has CUDA enabled but Omega_h does not...exiting")
  endif()
endif()

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


add_library(inclusion_solver STATIC inclusion_solver.cpp pfem_extras.cpp)
target_link_libraries(inclusion_solver mfem)
if(MFEM_USE_CUDA)
  enable_language(CUDA)
  set_source_files_properties(inclusion_solver.cpp pfem_extras.cpp PROPERTIES LANGUAGE CUDA)
  set_property(TARGET inclusion_solver PROPERTY CUDA_ARCHITECTURES ${INCLUSION_SOLVER_CUDA_ARCH})
endif()

add_executable(inclusion inclusion.cpp)
include_directories(${MFEM_INCLUDE_DIRS})
target_link_libraries(inclusion inclusion_solver)
if(MFEM_USE_CUDA)
  set_source_files_properties(inclusion.cpp PROPERTIES LANGUAGE CUDA)
  set_property(TARGET inclusion PROPERTY CUDA_ARCHITECTURES ${INCLUSION_SOLVER_CUDA_ARCH})
endif()

