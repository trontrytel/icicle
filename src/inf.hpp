/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref inf class - a facility for gathering simulation info
 */
#ifndef INF_HPP
#  define INF_HPP

#  include "cmn.hpp"
#  include "mtx.hpp"

#  ifdef USE_BOOST_TIMER
#    include <boost/timer/timer.hpp>
#  endif

#  include <sys/utsname.h> // uname 
#  include <unistd.h> // getlogin
#  include <map>

/// @brief a facility for gathering simulation info - timing, compiler configurration, etc
class inf : root
{
#  ifdef USE_BOOST_TIMER
  private: boost::timer::cpu_timer tmr;
#  endif
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
#  ifdef BZ_ARRAY_H
    im["blitz_version"] = string(BZ_VERSION);
    im["blitz_compiler_name"] = string(BZ__compiler_name);
    im["blitz_compiler_options"] = string(BZ__compiler_options);
    im["blitz_config_date"] = string(BZ__config_date);
    im["blitz_os_name"] = string(BZ__os_name);
    im["blitz_platform"] = string(BZ__platform);
#  endif

    // compile-time info defined in cfg.hpp
#  ifdef CONFIG_H
    im["compiler"] = string(INFO_COMPILER);
    im["compiler_version"] = string(INFO_COMPILER_VERSION);
    im["compiler_CXX_FLAGS"] = string(INFO_COMPILER_CXX_FLAGS);
    im["compiler_LINK_FLAGS"] = string(INFO_COMPILER_LINK_FLAGS);
    im["compiler_user"] = string(INFO_COMPILER_USER);
    im["compiler_system"] = string(INFO_COMPILER_SYSTEM);
    im["compiler_target"] = string(INFO_COMPILER_TARGET);
#  endif

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

#  ifdef USE_BOOST_TIMER
    {
      boost::timer::cpu_times t = tmr.elapsed(); // in nanoseconds
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
