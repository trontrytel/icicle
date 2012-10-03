set(tmp wkc_icmw8_case1_t0_test)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lboost_iostreams-mt -lboost_system-mt -lboost_filesystem-mt -lblitz")
add_test(${tmp} ${tmp})
