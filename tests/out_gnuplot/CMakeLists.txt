set(tmp test_out_gnuplot)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
add_test(${tmp} ${tmp})
