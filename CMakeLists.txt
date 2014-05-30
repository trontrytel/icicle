cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX C)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

# release-mode options
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Ofast -march=native") 

add_subdirectory(src) 

enable_testing()
add_subdirectory(tests)
