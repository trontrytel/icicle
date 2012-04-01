set(tmp test_todo_plot)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lplplotcxxd -lboost_iostreams-mt -lblitz")
add_test(${tmp} ${tmp})
