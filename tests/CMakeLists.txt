add_subdirectory(smolar_1983)
add_subdirectory(netcdf)

add_test(NAME courant_one_1D COMMAND ./courant_one_1D.sh)
