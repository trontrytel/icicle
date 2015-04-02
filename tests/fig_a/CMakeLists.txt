add_executable(calc calc.cpp)
add_test(calc calc ${CMAKE_BINARY_DIR})

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

foreach(micro blk_1m blk_2m)
  add_test(fig_a_diff_${micro} bash -c "
    for file in coord.h5 timestep0000000000.h5 timestep0000007200.h5 timestep0000009000.h5; do
      echo 'comparing $file ...'                                                               &&
      h5diff -v out_${micro}/$file ${CMAKE_CURRENT_SOURCE_DIR}/refdata/out_${micro}/$file || exit 1;
    done
  ")
endforeach()

foreach(micro blk_1m blk_2m lgrngn)
  add_executable(plot_${micro} plot_${micro}.cpp)
  add_test(plot_${micro} plot_${micro} ${CMAKE_BINARY_DIR})
  target_link_libraries(plot_${micro} ${HDF5_LIBRARIES})
  target_link_libraries(plot_${micro} ${Boost_LIBRARIES})
endforeach()

#add_executable(plot_lgrngn_spec plot_lgrngn_spec.cpp)
#add_test(plot_lgrngn_spec plot_lgrngn_spec ${CMAKE_BINARY_DIR})
#target_link_libraries(plot_lgrngn_spec ${HDF5_LIBRARIES})
#target_link_libraries(plot_lgrngn_spec ${Boost_LIBRARIES})
