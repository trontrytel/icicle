set(tmp test_todo_calc)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz -lcgicc -lboost_system-mt -lboost_filesystem-mt") # TODO! 

set(tmp test_todo_plot)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lboost_iostreams-mt -lboost_system-mt -lblitz") # TODO! 

set(tmp test_todo.cgi)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lcgicc -lboost_system -lboost_filesystem") # TODO! 
