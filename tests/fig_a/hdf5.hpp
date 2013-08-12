#pragma once

#include <blitz/array.h>
#include <H5Cpp.h>

blitz::Array<float, 2> h5load(
  const string &file, 
  const string &dataset
)
{
  notice_macro("about to open file: " << file)
  H5::H5File h5f(file, H5F_ACC_RDONLY);

  notice_macro("about to read dataset: " << dataset)
  H5::DataSet h5d = h5f.openDataSet(dataset);
  H5::DataSpace h5s = h5d.getSpace();

  if (h5s.getSimpleExtentNdims() !=4) 
    error_macro("need 3 dimensions and time")

  hsize_t n[4];
  enum {t, x, y, z};
  h5s.getSimpleExtentDims(n, NULL);

  blitz::Array<float, 2> tmp(n[x], n[z]);

  hsize_t 
    cnt[4] = { 1,      n[x], 1, n[z] }, 
    off[4] = { n[t]-1, 0,    0, 0    };
  h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);

  hsize_t ext[2] = {
    hsize_t(tmp.extent(0)), 
    hsize_t(tmp.extent(1))
  };
  h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(2, ext), h5s);

  return tmp;
}
