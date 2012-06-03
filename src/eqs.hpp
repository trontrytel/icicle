/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012 - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the @ref eqs class - a base class for all transport equation systems
 */
#ifndef EQS_HPP
#  define EQS_HPP

#  include "cmn.hpp"
#  include "rhs.hpp"

/// @brief a class defining a system of generalised transport equations
template <typename real_t>
class eqs 
{
  protected: struct groupid 
  {
    int id;
    groupid(int id=0) : id(id) {}
  };
  
  protected: struct nonnegative
  {
    bool is;
    nonnegative(bool is=false) : is(is) {}
  };

  /// @brief representation of a generalised transport equation 
  /// (e.g. eq. 19 in Smolarkiewicz & Margolin 1998)
  protected: struct gte {
    string name, desc, unit;
    nonnegative posdef;
    vector<int> pow_uvw; 
    groupid group;
    ptr_vector<rhs<real_t>> rhs_terms;
  };

  public: bool is_homogeneous(int e) 
  {
    assert(system().size() > e);
    if (system().at(e).rhs_terms.size() > 0) return false;
    return true;
  }

  protected: ptr_vector<struct eqs<real_t>::gte> sys;
  private: ptr_vector<gte> &system() { return sys; }

  public: int group(int e)
  {
    assert(system().size() > e);
    return system().at(e).group.id;
  }

  public: bool positive_definite(int e)
  {
    assert(system().size() > e);
    return system().at(e).posdef.is;
  }

  public: int n_vars()
  {
    return system().size();
  }

  public: int var_n_rhs_terms(int e)
  {
    assert(system().size() > e);
    return system().at(e).rhs_terms.size();
  }

  public: rhs<real_t> &var_rhs_term(int e, int i)
  {
    assert(system().size() > e);
    return system().at(e).rhs_terms[i];
  }

  public: bool var_dynamic(int e) // i.e. involved in calculation of velocities
  {
    assert(system().size() > e);
    if (system().at(e).pow_uvw.size() == 0) return false;
    assert(system().at(e).pow_uvw.size() == 3);
    return 
      system().at(e).pow_uvw[0] != 0 ||
      system().at(e).pow_uvw[1] != 0 ||
      system().at(e).pow_uvw[2] != 0;
  }
 
  public: string var_name(int i)
  {
    assert(system().size() > i);
    return system().at(i).name;
  }

  //       group   xyz   eq-pow      eq   pow       
  private: vector<vector<vector<pair<int, int>>>> vm;
  public: int n_group() { return vm.size(); }
  protected: void init_maps()
  {
    // mem allocation (vm) and group-related sanity checks
    {
      int n_groups = 1, group_last = -1;
      for (int e = 0; e < system().size(); ++e)
      {
        int group_curr = system().at(e).group.id;
        if (group_last != group_curr) 
        {
          assert(group_last + 1 == group_curr);
          n_groups++;
        }
        group_last = group_curr;
      }
      vm.resize(n_groups);
    }

    for (int g = 0; g < vm.size(); ++g)
    {
      vm[g].resize(3);
      for (int xyz = 0; xyz < 3; ++xyz)
      {
        // the first element is the one which will be multiplied/divided by others
        for (int e = 0; e < system().size(); ++e)
        {
          if (system().at(e).group.id != g) continue;
          assert(system().at(e).pow_uvw.size() == 3);
          if (system().at(e).pow_uvw[xyz] == 1)
          {
            assert(vm[g][xyz].size() == 0 && "TODO: some relevant message..."); // TODO: document this limitation...
            vm[g][xyz].push_back(pair<int,int>(e, 1));
          }
        }

        if (vm[g][xyz].size() == 0) continue;

        // and then goes the rest
        for (int e = 0; e < system().size(); ++e) 
        {
          if (system().at(e).group.id != g) continue;
          if (system().at(e).pow_uvw[xyz] != 0 && e != vm[g][xyz].begin()->first) 
            vm[g][xyz].push_back(pair<int,int>(e, system().at(e).pow_uvw[xyz]));
        }
      }
    }
  }

  public: vector<pair<int,int>> &velmap(int g, int xyz) // equation -> power
  {
    assert(vm.size() >= g);
    return vm[g][xyz];
  }
 
  public: string var_desc(int i)
  {
    assert(system().size() > i);
    return system().at(i).desc;
  }

  public: string var_unit(int i)
  {
    assert(system().size() > i);
    return system().at(i).unit;
  }

  protected: 
  template <class quan>
  string quan2str(quan q)
  {
    ostringstream tmp;
    tmp << boost::units::name_format << q;
    return tmp.str();
  }

  //////////////////////////
  /// auxiliary variables //
  //////////////////////////

  protected: struct invariable
  {
    bool is; 
    invariable(bool is=true) : is(is) {}
  };  

  public: struct axv { // TODO: should be private...
    string name, desc, unit; // TODO: name is already the key
    invariable constant;
    vector<int> dimspan; // TODO: document the meaning of the fourth dimenions
    vector<int> halo; // TODO: ditto
    static const int span = 0;
  };  
  
  protected: ptr_unordered_map<string, struct eqs<real_t>::axv> aux;
  public: ptr_unordered_map<string, struct eqs<real_t>::axv> &auxvars() { return aux; } // TODO: should be private!

  public: int n_auxv() { return auxvars().size(); } // TODO: is it needed anymore???

  public: mtx::idx aux_shape(const string &v, const mtx::idx &ijk) 
  {
    // sanity checks
    assert(auxvars()[v].dimspan.size() >= 3);
    assert(auxvars()[v].dimspan.size() >= auxvars()[v].halo.size());

    // computing the dimensions
    vector<mtx::rng> xyz(3);
    for (int d = 0; d < 3; ++d) 
    {
      int halo = auxvars()[v].halo.size() > d 
        ? auxvars()[v].halo[d] 
        : 0;
      xyz[d] = (auxvars()[v].dimspan[d] == axv::span) 
        ? mtx::rng(ijk[d].first() - halo, ijk[d].last() + halo)  
        : mtx::rng(0, auxvars()[v].dimspan[d] - 1);
    }

    return mtx::idx_ijk(xyz[0], xyz[1], xyz[2]);
  }

  public: vector<string> aux_names()
  {
    vector<string> names;
    for (const typename ptr_unordered_map<string, struct eqs<real_t>::axv>::value_type &axv : aux) names.push_back(axv.first);
    return names;
  }

  public: bool aux_const(const string &v)
  {
    assert(auxvars().count(v));
    return auxvars()[v].constant.is;
  } 

  public: bool aux_tobeoutput(const string &v)
  {
    assert(auxvars().count(v));
    return (!aux_const(v) 
      && auxvars()[v].dimspan[0] == axv::span
      && auxvars()[v].dimspan[0] == axv::span
      && auxvars()[v].dimspan[0] == axv::span);
  }

  public: string aux_desc(const string &v)
  {
    assert(auxvars().count(v));
    return auxvars()[v].desc;
  }

  public: string aux_unit(const string &v)
  {
    assert(auxvars().count(v));
    return auxvars()[v].unit;
  }

  /// @brief allows for post-advection and post-rhs adjustments to the state vector
  ///        (e.g. saturation adjustment in the ,,bulk'' cloud parameterisation
  ///        or some kind of smoothing)
  virtual void adjustments(
    int n,
    vector<ptr_vector<mtx::arr<real_t>>> &psi, // advected fields
    ptr_unordered_map<string, mtx::arr<real_t>> &aux, // auxiliary variables
    const ptr_vector<mtx::arr<real_t>> C, // Courant number fields
    const quantity<si::time, real_t> dt
  ) {} // no default adjustments

};
#endif
