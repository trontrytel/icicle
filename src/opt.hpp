/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - January 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_HPP
#  define OPT_HPP

#  include "cmn.hpp"

#  include <boost/program_options.hpp>
namespace po = boost::program_options;

template <typename real_t>
real_t real_cast(const po::variables_map& vm, const string &name) 
{
  string str;
  try 
  {
    str = vm[name].as<string>();
  }
  catch (boost::bad_any_cast &)
    error_macro("option " << name << " not specified")

  try
  {
    return boost::lexical_cast<real_t>(vm[name].as<string>());
  }
  catch (boost::bad_lexical_cast &) 
    error_macro("failed to convert " << vm[name].as<string>() 
      << " into floating point value while handling the --" << name << " option"
  )
}

#endif
