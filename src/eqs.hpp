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
    string name;
    string desc;
    string unit; 
    //vector<rhs> // source terms - TODO
  };

  public: virtual vector<gte> &system() = 0;
  public: virtual bool has_qx() { return false; }
  public: virtual bool has_qy() { return false; }
  public: virtual bool has_qz() { return false; }
  public: virtual int idx_qx() { assert(false); }
  public: virtual int idx_qy() { assert(false); }
  public: virtual int idx_qz() { assert(false); }

  public: int n_vars()
  {
    return system().size();
  }
 
  public: string var_name(int i)
  {
    // TODO try/catch if i within range
    return system().at(i).name;
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
