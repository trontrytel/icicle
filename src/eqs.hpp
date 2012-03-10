/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the @ref eqs class - a base class for all transport equation systems
 */
#ifndef EQS_HPP
#  define EQS_HPP

#  include "cmn.hpp"
#  include "rhs.hpp"

/// @brief a container for systems of generalised transport equations
template <typename real_t>
class eqs : root
{
  // TODO: ctor with sanity checks for name uniqueness?

  /// @brief representation of a generalised transport equation 
  /// (e.g. eq. 19 in Smolarkiewicz & Margolin 1998)
  protected: class groupid {public:int id;public:groupid(int id=0):id(id){}}; // yeah!
  protected: class positive_definite {public:bool is;public:positive_definite(bool is=false):is(is){}}; // yeah!
  protected: struct gte {
    string name, desc, unit;
    positive_definite posdef;
    vector<int> pow_uvw; 
    groupid group;
    ptr_vector<rhs<real_t>> rhs_terms;
  };

  public: bool is_homogeneous(int e) {
    if (system().at(e).rhs_terms.size() > 0) return false;
    return true;
  }

  public: virtual ptr_vector<gte> &system() = 0;

  public: int group(int e)
  {
    return system().at(e).group.id;
  }

  public: bool positive_definite(int e)
  {
    return system().at(e).posdef.is;
  }

  public: int n_vars()
  {
    return system().size();
  }

  public: int var_n_rhs_terms(int e)
  {
    return system().at(e).rhs_terms.size();
  }

  public: rhs<real_t> &var_rhs_term(int e, int i)
  {
    return system().at(e).rhs_terms[i];
  }

  public: bool var_dynamic(int e) // i.e. involved in calculation of velocities
  {
    if (system().at(e).pow_uvw.size() == 0) return false;
    assert(system().at(e).pow_uvw.size() == 3);
    return 
      system().at(e).pow_uvw[0] != 0 ||
      system().at(e).pow_uvw[1] != 0 ||
      system().at(e).pow_uvw[2] != 0;
  }
 
  public: string var_name(int i)
  {
    // TODO try/catch if i within range
    return system().at(i).name;
  }

  //       group   xyz   eq-pow      eq   pow       
  private: vector<vector<vector<pair<int, int>>>> vm;
  private: void init_maps()
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
    if (vm.size() == 0) init_maps();
    return vm[g][xyz];
  }
 
  public: string var_desc(int i)
  {
    // TODO try/catch if i within range
    return system().at(i).desc;
  }

  public: string var_unit(int i)
  {
    // TODO try/catch if i within range
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
};
#endif
