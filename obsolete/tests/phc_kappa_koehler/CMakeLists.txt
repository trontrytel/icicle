set(tmp phc_kappa_koehler_test)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
add_test(${tmp} ${tmp})
