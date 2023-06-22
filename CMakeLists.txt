cmake_minimum_required(VERSION 3.5)

project("FENGSHA CCPP"
        VERSION 0.1.0
        LANGUAGES Fortran)

add_subdirectory(src)

add_subdirectory(driver/ccpp)
