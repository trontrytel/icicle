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

template <typename real_t>
class eqs : root
{
  // TODO: ctor with sanity checks for name uniqueness?

  // A generalised transport equation (e.g. eq. 19 in Smolarkiewicz & Margolin 1998)
  protected: struct gte {
    string name, desc, unit;
    int pow_u, pow_v, pow_w;
    //vector<rhs> // source terms - TODO
  };

  public: virtual vector<gte> &system() = 0;

  public: int n_vars()
  {
    return system().size();
  }

  public: bool var_dynamic(int i) // i.e. involved in calculation of velocities
  {
    return 
      system().at(i).pow_u != 0 ||
      system().at(i).pow_v != 0 ||
      system().at(i).pow_w != 0;
  }
 
  public: string var_name(int i)
  {
    // TODO try/catch if i within range
    return system().at(i).name;
  }
 
  // TODO shirten it!
  public: map<int,int> velmap_x() // equeation -> power
  {
    map<int, int> m;
    for (int e = 0; e < system().size(); ++e)
    {
      if (system().at(e).pow_u != 0)
      {
        m[e] = system().at(e).pow_u;
        cerr << "mx[" << e << "]=" << m[e] << endl;
      }
    }
    return m;
  }
  public: map<int,int> velmap_y() // equeation -> power
  {
    map<int, int> m;
    for (int e = 0; e < system().size(); ++e)
      if (system().at(e).pow_v != 0)
      {
        m[e] = system().at(e).pow_v;
        cerr << "my[" << e << "]=" << m[e] << endl;
      }
    return m;
  }
  public: map<int,int> velmap_z() // equeation -> power
  {
    map<int, int> m;
    for (int e = 0; e < system().size(); ++e)
      if (system().at(e).pow_w != 0) m[e] = system().at(e).pow_w;
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
