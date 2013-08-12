#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "common.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string opts_common = 
    "--outfreq=10 --nt 10";
  set<string> opts_micro({
    "--micro=blk_1m --outfile=out_blk_1m.h5",
    "--micro=blk_2m --outfile=out_blk_2m.h5",
    "--micro=lgrngn --outfile=out_lgrngn.h5 --backend=CUDA --sd_conc_mean=64 --sstp_cond=20 --RH_max 1.2" // TODO!
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/icicle " << opts_common << " " << opts_m;  

    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}
