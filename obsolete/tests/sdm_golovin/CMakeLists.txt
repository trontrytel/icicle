find_package(Boost COMPONENTS system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(tmp sdm_golovin_test)
add_executable(${tmp} ${tmp}.cpp)
add_dependencies(${tmp} icicle)
target_link_libraries(${tmp} ${Boost_LIBRARIES})
add_test(${tmp} ${tmp})
