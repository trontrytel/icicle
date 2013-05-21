find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(tmp calc_chem)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz -lcgicc") # TODO! 
target_link_libraries(${tmp} ${Boost_LIBRARIES})

set(tmp plot_chem)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz") # TODO! 
target_link_libraries(${tmp} ${Boost_LIBRARIES})

#set(tmp plot_spec)
#add_executable(${tmp} ${tmp}.cpp)
#add_dependencies(${tmp} icicle)
#target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz") # TODO! 
#target_link_libraries(${tmp} ${Boost_LIBRARIES})
