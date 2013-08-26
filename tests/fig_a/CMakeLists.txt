add_executable(calc calc.cpp)
add_test(calc calc ${CMAKE_BINARY_DIR})

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

foreach(micro blk_1m blk_2m lgrngn)
  add_executable(plot_${micro} plot_${micro}.cpp)
  add_test(plot_${micro} plot_${micro} ${CMAKE_BINARY_DIR})
  target_link_libraries(plot_${micro} ${HDF5_LIBRARIES})
  target_link_libraries(plot_${micro} ${Boost_LIBRARIES})
endforeach()

add_executable(plot_lgrngn_spec plot_lgrngn_spec.cpp)
add_test(plot_lgrngn_spec plot_lgrngn_spec ${CMAKE_BINARY_DIR})
target_link_libraries(plot_lgrngn_spec ${HDF5_LIBRARIES})
target_link_libraries(plot_lgrngn_spec ${Boost_LIBRARIES})

