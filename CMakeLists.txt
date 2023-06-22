cmake_minimum_required(VERSION 3.5)

project("FENGSHA CCPP"
    VERSION 0.1.0
    LANGUAGES Fortran)

# If GNU compiler and Debug build type, add extra debug flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95" AND CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-g -Wall -Wextra -Og -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -fmax-errors=5)
endif()

add_subdirectory(src)

get_target_property(CFLAGS "FENGSHA" COMPILE_OPTIONS)
if (CFLAGS AND NOT "${CFLAGS}" STREQUAL "")
    message(STATUS "FENGHSA Fortran compiler flags: ${CFLAGS}")
endif()

add_subdirectory(driver/ccpp)
