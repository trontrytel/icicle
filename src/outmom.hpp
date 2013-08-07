#pragma once

template <typename real_t>
using outmom_t = std::map<std::pair<
  quantity<si::length, real_t>,
  quantity<si::length, real_t>
>, std::vector<int> >;

// TODO: option parsing code should go here...
