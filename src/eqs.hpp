/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_HPP
#  define EQS_HPP

#  include "cmn.hpp" // root class, error reporting

template <typename real_t>
class eqs : root
{
  // A generalised transport equation (e.g. eq. 19 in Smolarkiewicz & Margolin 1998)
  protected: struct gte {
    string name;
    string unit; 
    //vector<rhs> // source terms - TODO
  };

  public: virtual vector<gte> &system() = 0;
  // TODO: jacobian() ?
  // TODO: particles() ?

  public: int n_vars()
  {
    return system().size();
  }
 
  public: string var_name(int i)
  {
    // TODO try/catch if i within range
    return system().at(i).name;
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
