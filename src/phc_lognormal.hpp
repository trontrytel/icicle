/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date May 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef PHC_LOGNORMAL_HPP
#  define PHC_LOGNORMAL_HPP

namespace phc
{
  // lognormal distribution as a function of ln(r) (Seinfeld & Pandis 1997 eq 7.33)
  phc_declare_funct_macro quantity<power_typeof_helper<si::length,static_rational<-3>>::type, real_t> log_norm_n_e(
    quantity<si::length, real_t> mean_r,
    quantity<si::dimensionless, real_t> stdev, 
    quantity<power_typeof_helper<si::length,static_rational<-3>>::type, real_t> n_tot, 
    quantity<si::dimensionless, real_t> lnr
  )
  {
    return n_tot / sqrt(2*M_PI) / log(stdev) * exp(-pow((lnr - log(mean_r/si::metres)), 2) / 2 / pow(log(stdev),2));
  }

  // lognormal distribution as a function of r (Seinfeld & Pandis 1997 eq 7.34)
  phc_declare_funct_macro quantity<power_typeof_helper<si::length,static_rational<-4>>::type, real_t> log_norm_n(
    quantity<si::length, real_t> mean_r,
    quantity<si::dimensionless, real_t> stdev, 
    quantity<power_typeof_helper<si::length,static_rational<-3>>::type, real_t> n_tot, 
    quantity<si::length, real_t> r
  )
  {
    return n_tot / sqrt(2*M_PI) / log(stdev) * exp(-pow((log(r/mean_r)), 2) / 2 / pow(log(stdev),2)) / r;
  }
};

#endif
