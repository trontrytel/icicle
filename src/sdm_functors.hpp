/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief a set of functors used with Thrust routines in the eqs_todo_sdm class
 */
#pragma once

#include "sdm_base.hpp"
#include "phc_lognormal.hpp"

#ifdef USE_THRUST

namespace sdm {

  /// @brief a random-number-generator functor
  // TODO: random seed as an option
  // TODO: RNG engine as an option
  template <typename real_t> 
  class rng 
  {
    private: thrust::random::taus88 engine; 
    private: thrust::uniform_real_distribution<real_t> dist;
    public: rng(real_t a, real_t b, real_t seed) : dist(a, b), engine(seed) {}
    public: real_t operator()() { return dist(engine); }
  };

  /// @brief a functor that divides by real constant and cast to int
  template <typename real_t> 
  class divide_by_constant
  {
    private: real_t c;
    public: divide_by_constant(real_t c) : c(c) {}
    public: int operator()(real_t x) { return x/c; }
  };

  /// @brief a functor that multiplies by a real constant
  template <typename real_t> 
  class multiply_by_constant
  {
    private: real_t c;
    public: multiply_by_constant(real_t c) : c(c) {}
    public: real_t operator()(real_t x) { return x*c; }
  };

  /// @brief a functor that ravels (i,j) index pairs into a single index
  class ravel_indices
  {
    private: int n;
    public: ravel_indices(int n) : n(n) {}
    public: int operator()(int i, int j) { return i + j * n; }
  };

  /// @brief a functor interface to fmod()
  template <typename real_t> 
  class modulo
  {
    private: real_t mod;
    public: modulo(real_t mod) : mod(mod) {}
    public: real_t operator()(real_t a) { return fmod(a + mod, mod); }
  };

  /// @brief a Thrust-to-Blitz data transfer functor 
  template <typename real_t> 
  class thrust2blitz
  {
    private: int n;
    private: const thrust::device_vector<int> &idx2ij;
    private: const thrust::device_vector<real_t> &from;
    private: mtx::arr<real_t> &to;
    private: real_t scl;

    // ctor
    public: thrust2blitz(int n, 
      const thrust::device_vector<int> &idx2ij,
      const thrust::device_vector<real_t> &from,
      mtx::arr<real_t> &to,
      real_t scl = real_t(1)
    ) : n(n), idx2ij(idx2ij), from(from), to(to), scl(scl) {}

    public: void operator()(int idx) 
    { 
      to(idx2ij[idx] % n, idx2ij[idx] / n, 0) = scl * from[idx]; 
    }
  };

  /// @brief a Blitz-to-Thrust data transfer functor
  template <typename blitz_real_t, typename real_t>  // TODO: blitz_real_t == real_t ? 
  class blitz2thrust
  {
    private: int n;
    private: const mtx::arr<blitz_real_t> &from;
    private: thrust::device_vector<real_t> &to;
    private: real_t scl;

    // ctor
    public: blitz2thrust(int n,
      const mtx::arr<blitz_real_t> &from,
      thrust::device_vector<real_t> &to,
      real_t scl = real_t(1)
    ) : n(n), from(from), to(to), scl(scl) {}

    public: void operator()(int ij)
    {
      to[ij] = scl * from(ij % n, ij / n, 0);
    }
  };

  /// @brief a functor interface to phc_lognormal.hpp routines
  template <typename real_t> 
  class lognormal
  {
    private: quantity<si::length, real_t> mean_rd1, mean_rd2;
    private: quantity<si::dimensionless, real_t> sdev_rd1, sdev_rd2;
    private: quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_tot, n2_tot;
    
    public: lognormal(
      const quantity<si::length, real_t> mean_rd1,
      const quantity<si::dimensionless, real_t> sdev_rd1,
      const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_tot,
      const quantity<si::length, real_t> mean_rd2,
      const quantity<si::dimensionless, real_t> sdev_rd2,
      const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n2_tot
    ) : mean_rd1(mean_rd1), sdev_rd1(sdev_rd1), n1_tot(n1_tot), mean_rd2(mean_rd2), sdev_rd2(sdev_rd2), n2_tot(n2_tot) {}

    public: real_t operator()(real_t lnrd)
    {
      return real_t(( //TODO allow more modes of distribution
// TODO: logarith or not: as an option
          phc::log_norm_n_e<real_t>(mean_rd1, sdev_rd1, n1_tot, lnrd) + 
          phc::log_norm_n_e<real_t>(mean_rd2, sdev_rd2, n2_tot, lnrd) 
        ) * si::cubic_metres
      );
    }
  };

  template <typename real_t, class xi>
  class moment_counter
  { 
    private: const sdm::stat_t<real_t> &stat;
    private: const real_t threshold;
    private: int k;
    public: moment_counter(
      const sdm::stat_t<real_t> &stat,
      const real_t threshold,
      const int k
    ) : stat(stat), threshold(xi::xi_of_rw(threshold)), k(k) {}
    public: real_t operator()(const thrust_size_t id) const
    {   
      switch (k)
      {
        case 0: return stat.xi[id] > threshold ? (stat.n[id]) : 0;
        case 1: return stat.xi[id] > threshold ? (stat.n[id] * xi::rw_of_xi(stat.xi[id])) : 0; 
        case 2: return stat.xi[id] > threshold ? (stat.n[id] * xi::rw2_of_xi(stat.xi[id])) : 0;
        case 3: return stat.xi[id] > threshold ? (stat.n[id] * xi::rw3_of_xi(stat.xi[id])) : 0;
        default: return stat.xi[id] > threshold ? (stat.n[id] * pow(xi::rw_of_xi(stat.xi[id]), k)) : 0;
      }
    }
  };

  template <typename real_t>
  class chem_counter 
  { 
    private: const sdm::stat_t<real_t> &stat;
    private: const int spec;
    public: chem_counter(
      const sdm::stat_t<real_t> &stat,
      const int spec
    ) : stat(stat), spec(spec) {}
    public: real_t operator()(const thrust_size_t id) const
    {   
      return stat.n[id] * stat.c_aq[spec * stat.n_part + id]; // TODO: assumed c_aq
    }
  };

};
#endif
