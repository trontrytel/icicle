add_executable(calc_b calc.cpp)
add_test(calc_b calc ${CMAKE_BINARY_DIR})

find_package(Boost COMPONENTS system timer REQUIRED)
target_link_libraries(calc_b ${Boost_LIBRARIES})

#
#foreach(micro blk_1m blk_2m lgrngn)
#  add_executable(plot_${micro} plot_${micro}.cpp)
#  add_test(plot_${micro} plot_${micro} ${CMAKE_BINARY_DIR})
#  target_link_libraries(plot_${micro} ${HDF5_LIBRARIES})
#endforeach()
#
#add_executable(plot_lgrngn_spec plot_lgrngn_spec.cpp)
#add_test(plot_lgrngn_spec plot_lgrngn_spec ${CMAKE_BINARY_DIR})
#target_link_libraries(plot_lgrngn_spec ${HDF5_LIBRARIES})
#target_link_libraries(plot_lgrngn_spec ${Boost_LIBRARIES})
