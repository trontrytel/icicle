#pragma once
// Zach's mass-doubling bin layout

vector<quantity<si::length, real_t>> bins()
{
  const int n_bins = 48;
  vector<quantity<si::length, real_t>> bins(n_bins);
  quantity<si::length, real_t> r = real_t(.5e-8) * si::metres;
  for (int i = 0; i < n_bins; ++i)
  {
    bins[i] = r;
    quantity<si::volume, real_t> v = real_t(4./3) * phc::pi<real_t>() * pow<3>(r);
    r = root<3>(real_t(3./4) / phc::pi<real_t>() * 2 * v);
  }
  return bins;
}
