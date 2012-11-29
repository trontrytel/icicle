find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(tmp wkc_icmw8_case1_t0_test)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} "-lnetcdf_c++4 -lblitz")
target_link_libraries(${tmp} ${Boost_LIBRARIES})
add_dependencies(${tmp} icicle)
add_test(${tmp} ${tmp})
