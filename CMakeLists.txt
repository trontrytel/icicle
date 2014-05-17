cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX C)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

add_subdirectory(src) 

enable_testing()
add_subdirectory(tests)
