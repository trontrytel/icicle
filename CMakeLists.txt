cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX C)

add_subdirectory(src) 

# TODO: automate "make doc", Some hints here: http://www.cmake.org/pipermail/cmake/2007-February/012803.html

# tests
enable_testing()
add_subdirectory(tests)
