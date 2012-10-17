set(tmp sdm_golovin_test)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
add_test(${tmp} ${tmp})
