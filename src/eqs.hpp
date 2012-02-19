/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_HPP
#  define EQS_HPP

#  include "cmn.hpp"
#  include "rhs.hpp"

template <typename real_t>
class eqs : root
{
  // TODO: ctor with sanity checks for name uniqueness?

  // A generalised transport equation (e.g. eq. 19 in Smolarkiewicz & Margolin 1998)
  protected: struct gte {
    string name, desc, unit;
    vector<int> pow_uvw;
    ptr_vector<rhs<real_t>> source_terms;
  };

  public: bool is_homogeneous() {
    for (int e = 0; e < system().size(); ++e) 
      if (system().at(e).source_terms.size() > 0) return false;
    return true;
  }

  public: virtual ptr_vector<gte> &system() = 0;

  public: int n_vars()
  {
    return system().size();
  }

  public: int var_n_rhs(int e)
  {
    return system().at(e).source_terms.size();
  }

  public: rhs<real_t> &var_rhs(int e, int i)
  {
    return system().at(e).source_terms[i];
  }

  public: bool var_dynamic(int e) // i.e. involved in calculation of velocities
  {
    if (system().at(e).pow_uvw.size() == 0) return false;
    assert(system().at(e).pow_uvw.size() == system().size());
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
 
  public: map<int,int> velmap(int xyz) // equation -> power
  {
    map<int, int> m;
    for (int e = 0; e < system().size(); ++e)
    {
      assert(system().at(e).pow_uvw.size() == system().size());
      if (system().at(e).pow_uvw[xyz] != 0) 
        m[e] = system().at(e).pow_uvw[xyz];
    }
    return m;
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
