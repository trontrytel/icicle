set(tmp adv_courant_one_test)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
add_test(${tmp} ${tmp})
