cmake_minimum_required (VERSION 2.8)
project (SHARPSAT) 

option(STATIC "STATIC" OFF)
set(CMAKE_VERBOSE_MAKEFILE on)

if(STATIC)
set(CMAKE_EXE_LINKER_FLAGS -static)
set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
set(CMAKE_SHARED_LIBRARY_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
set(BUILD_SHARED_LIBS OFF)
endif(STATIC)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra -Wshadow -g -Wfatal-errors -Wno-sign-compare -Wno-literal-suffix")

SET(CMAKE_C_FLAGS "-std=c99 -Wall -msse4.1 -mpclmul")

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)

add_subdirectory(src)
add_subdirectory(flowcutter)

add_executable (sharpSATU src/main.cpp)
if(STATIC)
SET_TARGET_PROPERTIES (sharpSATU PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)
target_link_libraries (sharpSATU ${GMPXX_LIB} ${GMP_LIB} libsharpSAT)

add_executable (sharpSATW src/main.cpp)
if(STATIC)
SET_TARGET_PROPERTIES (sharpSATW PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)
SET_TARGET_PROPERTIES(sharpSATW PROPERTIES COMPILE_FLAGS "-DWEIGHTED")
target_link_libraries (sharpSATW ${GMPXX_LIB} ${GMP_LIB} libsharpSAT)

add_executable (flow_cutter_pace17 flowcutter/pace.cpp)
if(STATIC)
SET_TARGET_PROPERTIES (flow_cutter_pace17 PROPERTIES LINK_SEARCH_START_STATIC 1)
endif(STATIC)
target_link_libraries (flow_cutter_pace17 libflowcutter)