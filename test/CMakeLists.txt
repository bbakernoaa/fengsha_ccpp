add_executable(FENGSHA_test
    test.f90
)

target_include_directories(FENGSHA_test
    PRIVATE $<TARGET_PROPERTY:FENGSHA,BINARY_DIR>
)
target_link_libraries(FENGSHA_test
    PRIVATE FENGSHA
)
add_test(NAME Basic COMMAND FENGSHA_test)
