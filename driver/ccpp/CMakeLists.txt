add_library(FENGSHA_CCPP_WRAPPER
    catchem_dust_wrapper.F90
)

target_link_libraries(FENGSHA_CCPP_WRAPPER
    PUBLIC FENGSHA
)
target_include_directories(FENGSHA_CCPP_WRAPPER
    PUBLIC $<TARGET_PROPERTY:FENGSHA,BINARY_DIR>
)
