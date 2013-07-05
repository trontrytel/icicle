
vector<quantity<si::length, real_t>> bins_dry()
{
  vector<quantity<si::length, real_t>> ret;
  // dry radius bins: .001  ...  .01  ...  .1  ... 1 (30 bins in total)
  for (int i = 0; i < 30; ++i) 
    ret.push_back(real_t(1e-6 * pow(10, -3 +   i   * .1)) * si::metres);
  return ret;
}

vector<quantity<si::length, real_t>> bins_wet()
{
  vector<quantity<si::length, real_t>> ret;
  for (int i = 0; i < 50; ++i) 
    ret.push_back(real_t(1e-6 * pow(10, -3 +   i   * .1)) * si::metres);
  return ret;
}
