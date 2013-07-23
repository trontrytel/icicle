/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref inf class - a facility for gathering simulation info
 */
#include "inf.hpp"

#include <sys/utsname.h> // uname 
#include <unistd.h> // getlogin

#include "cmn/cmn_error.hpp"

#include "cfg/cfg_info.hpp"

#include <sstream>
using std::ostringstream;

map<string,string> inf::get_map()
{
  map<string,string> im;
  im["command_line"] = options;
  im["built_date"] = string(__DATE__);

  // Blitz++ settings
  im["blitz_version"] = string(BZ_VERSION);
  im["blitz_compiler_name"] = string(BZ__compiler_name);
  im["blitz_compiler_options"] = string(BZ__compiler_options);
  im["blitz_config_date"] = string(BZ__config_date);
  im["blitz_os_name"] = string(BZ__os_name);
  im["blitz_platform"] = string(BZ__platform);

  // compile-time info defined in cfg.hpp
  im["compiler"] = string(INFO_COMPILER);
  im["compiler_version"] = string(INFO_COMPILER_VERSION);
  im["compiler_CXX_FLAGS"] = string(INFO_COMPILER_CXX_FLAGS);
  im["compiler_LINK_FLAGS"] = string(INFO_COMPILER_LINK_FLAGS);
  im["compiler_user"] = string(INFO_COMPILER_USER);
  im["compiler_system"] = string(INFO_COMPILER_SYSTEM);
  im["compiler_target"] = string(INFO_COMPILER_TARGET);

  // run-time info 
  {   
    struct utsname uts;
    uname(&uts);
    im["runtime_system"] = string(uts.sysname) + " " + string(uts.machine);
    char* login = getlogin();
    if (login == NULL) warning_macro("getlogin() failed!") // happens e.g. with the "screen" terminal
    string tmp = (login == NULL ? "?" : login) + string("@") + uts.nodename;
    im["runtime_user"] =tmp;
  }   

  return im;
}
