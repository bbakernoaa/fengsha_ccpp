set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(FENGSHA SHARED
    catchem_constants.F90
    catchem_config.F90
    dust_data_mod.F90
    dust_fengsha_mod.F90
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_BINARY_DIR})
