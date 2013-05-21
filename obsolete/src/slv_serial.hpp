/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "slv.hpp"
#include "stp.hpp"
#include "tmp.hpp" // TODO: get rid of it!

#  include <memory>
using std::unique_ptr;

template <typename real_t>
class slv_serial : public slv<real_t>
{
  private: ptr_vector<typename adv<real_t>::op3D> advop; // advop[positive_definite]
  private: const mtx::idx_ijk ijk;
  private: out<real_t> &output;
  private: const stp<real_t> &setup;
  private: vector<int> halo_sclr; // TODO: within ctor?
  private: int halo_vctr; // TODO: within ctor?

  // Courant numbers (Cx,Cy,Cz)
  private: ptr_vector<mtx::arr<real_t>> C; // TODO: nullable? (uwaga na kod w FCT! ii = i+/-1

  // advected vector fields at n time levels (2 for MPDATA, 3 for leapfrog, etc)
  private: vector<ptr_vector<mtx::arr<real_t>>> psi;

  // auxiliary variables 
  private: ptr_unordered_map<string, mtx::arr<real_t>> aux;

  // RHS representations (implicit and explicit parts)
  private: ptr_vector<boost::nullable<mtx::arr<real_t>>> rhs_R, rhs_C;

  private: unique_ptr<tmp<real_t>> cache; 
  private: vector<mtx::arr<real_t>*> tmpvec; // used in update_forcings
  private: vector<ptr_vector<mtx::arr<real_t>>> Q; // used in update_courants

  // ctor
  public: slv_serial(
    const stp<real_t> &setup, 
    out<real_t> &output,
    int i_min, int i_max,
    int j_min, int j_max,
    int k_min, int k_max
  );

  public: void copy(int from, int to);

  public: void fill_with_nans(int e, int n);

  public: void record(const int n, const unsigned long t);

  public: typename mtx::arr<real_t>::type data(int e, int n, const mtx::idx &idx);

  public: void fill_halos(const int e, const int n);

  private:
  template<class idx>
  void fill_halos_helper(const int nghbr, const int e, const int n, 
    const int i_min, const int i_max, const mtx::rng &j, const mtx::rng &k, int mod
  );

  public: void advect(int e, int n, int s);

  public: void update_courants(const int g, const int nm1, const int nm0);

  /// \param n the time level to use for updating the forcings
  public: void update_forcings(int n /*, const quantity<si::time, real_t> t*/);

  public: void apply_forcings(int e, int n, quantity<si::time, real_t> dt);

  public: void apply_adjustments(int n, const quantity<si::time, real_t> dt, bool record);

  public: void cycle_arrays(const int e, const int n); // TODO: n unneeded?

  private: mtx::arr<real_t> *stash = NULL;
  private: bool stash_empty = true;
  public: void stash_cycle(int e, int n);

  public: void integ_loop();
};
