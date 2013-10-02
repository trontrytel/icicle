#pragma once

// bin sizes for calc and plot
vector<quantity<si::length>> bins_dry()
{
  vector<quantity<si::length>> ret;
  // dry radius bins: .001 ... .01 ... .1 ... 1 (30 bins in total)
  for (int i = 0; i < 30; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .1) * si::metres);
  return ret;
}

vector<quantity<si::length>> bins_wet()
{
  vector<quantity<si::length>> ret;
  for (int i = 0; i < 50; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .2) * si::metres);
  return ret;
}

// focus plot locations
int ox = 0, oy=1;
std::pair<
  std::set<std::pair<int,int>>,
  std::set<std::pair<int,int>>
> focus = {
  {   // left column
    {17+ox, 24+oy}, 
    {17+ox, 35+oy}, 
    {17+ox, 37+oy}, 
    {17+ox, 41+oy},
    {17+ox, 70+oy} 
  },{ // right column
    {58+ox, 24+oy}, 
    {58+ox, 35+oy}, 
    {58+ox, 37+oy}, 
    {58+ox, 41+oy},
    {58+ox, 70+oy}
  }
};
