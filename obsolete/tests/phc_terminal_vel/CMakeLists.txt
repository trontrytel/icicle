find_package(Boost COMPONENTS iostreams system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(tmp phc_terminal_vel_test)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz")
target_link_libraries(${tmp} ${Boost_LIBRARIES})
add_dependencies(${tmp} icicle)
add_test(${tmp} ${tmp})
