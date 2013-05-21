cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX C)

add_subdirectory(src) 

# TODO: automate "make doc"

# tests
enable_testing()
add_subdirectory(tests)
