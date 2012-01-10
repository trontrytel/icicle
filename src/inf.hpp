/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INF_HPP
#  define INF_HPP

#  include "cmn.hpp"

#  ifdef USE_BOOST_TIMER
#    include <boost/timer/timer.hpp>
using namespace boost::timer;
#  endif

#  include <sys/utsname.h> // uname 
#  include <unistd.h> // getlogin
#  include <map>

class inf : root
{
  private: cpu_timer tmr;
  private: string options;
  public: inf(const string &options)
    : options(options)
  { }

  public: map<string,string> get_map()
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
      if (login == NULL) warning_macro("getlogin() failed!")
      string tmp = (login == NULL ? "?" : login) + string("@") + uts.nodename;
      im["runtime_user"] =tmp;
    }   

#  ifdef USE_BOOST_TIMER
    {
      cpu_times t = tmr.elapsed(); // in nanoseconds
      {
        ostringstream tmp;
        tmp << double(t.wall) * 1e-9;
        im["wall_time_seconds"] = tmp.str();
      }
      {
        ostringstream tmp;
        tmp << double(t.user) * 1e-9;
        im["user_time_seconds"] = tmp.str();
      }
      {
        ostringstream tmp;
        tmp << double(t.system) * 1e-9;
        im["system_time_seconds"] = tmp.str();
      }
    }
#  endif
    return im;
  }
};

#endif
