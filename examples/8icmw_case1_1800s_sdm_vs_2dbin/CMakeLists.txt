find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(tmp calc)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz") # TODO! 
target_link_libraries(${tmp} ${Boost_LIBRARIES})

set(tmp plot)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz") # TODO! 
target_link_libraries(${tmp} ${Boost_LIBRARIES})
