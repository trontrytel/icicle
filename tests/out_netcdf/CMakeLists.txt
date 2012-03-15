set(tmp test_out_netcdf)
add_executable(${tmp} ${tmp}.cpp)
add_test(${tmp} ${tmp})
