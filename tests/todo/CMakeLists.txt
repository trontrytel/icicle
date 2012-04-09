set(tmp test_todo_calc)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz") # TODO! 
add_test(${tmp} ${tmp})

set(tmp test_todo_plot)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lboost_iostreams-mt -lboost_system-mt -lblitz") # TODO! 
add_test(${tmp} ${tmp})
