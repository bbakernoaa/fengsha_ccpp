add_library(FENGSHA_CCPP_WRAPPER
    catchem_dust_wrapper.F90
)

target_include_directories(FENGSHA_CCPP_WRAPPER
    PRIVATE $<TARGET_PROPERTY:FENGSHA,BINARY_DIR>
)
