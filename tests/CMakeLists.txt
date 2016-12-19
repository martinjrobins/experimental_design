
set(test_package 
    e_surface 
    )

include_directories(${sinusoidal_voltammetry_SOURCE_DIR}/src)

foreach(test ${test_package})
    cxxtest_add_test(${test} ${test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${test}.hpp)
    target_link_libraries(${test} sinusoidal_voltammetry)
endforeach(test)
