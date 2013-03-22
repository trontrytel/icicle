/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "opt.hpp"
#include "eqs/eqs_scalar_advection.hpp"
#include "eqs/eqs_isentropic.hpp"
#include "eqs/eqs_todo_bulk_ode.hpp"
#include "eqs/eqs_todo_mm.hpp"
#include "eqs/eqs_todo_sdm.hpp"

#include <boost/spirit/include/qi.hpp>        
#include <boost/fusion/adapted/std_pair.hpp>  

#include <list>

inline void opt_eqs_desc(po::options_description &desc)
{
  desc.add_options()
    ("eqs", po::value<string>(), "equation system: todo, isentropic, ...")

    ("eqs.isentropic.nlev", po::value<int>(), "number of fluid layers")
    ("eqs.isentropic.p_top", po::value<string>()->default_value("0"), "pressure at the uppermost surface [Pa]")
    ("eqs.isentropic.theta_frst", po::value<string>(), "mid-layer potential temperature of the first layer [K]")
    ("eqs.isentropic.abslev", po::value<int>(), "absorber lowermost level")
    ("eqs.isentropic.absamp", po::value<string>()->default_value("1"), "absorber amplitude [1]") 

    ("eqs.todo_bulk.cond", po::value<bool>(), "cloud water condensation [on/off]")
    ("eqs.todo_bulk.cevp", po::value<bool>(), "cloud water evaporation [on/off]")
    ("eqs.todo_bulk.conv", po::value<bool>(), "conversion of cloud water into rain [on/off]")
    ("eqs.todo_bulk.clct", po::value<bool>(), "collection of cloud water by rain [on/off]")
    ("eqs.todo_bulk.sedi", po::value<bool>(), "rain water sedimentation [on/off]")
    ("eqs.todo_bulk.revp", po::value<bool>(), "rain water evaporation [on/off]")

    ("eqs.todo_mm.act",  po::value<bool>(), "cloud water activation [on/off]")
    ("eqs.todo_mm.cond", po::value<bool>(), "condensation/evaporation [on/off]")
    ("eqs.todo_mm.acc",  po::value<bool>(), "collection of cloud water by rain [on/off]")
    ("eqs.todo_mm.autoc",po::value<bool>(), "autoconversion of cloud water by rain [on/off]")
    ("eqs.todo_mm.turb", po::value<bool>(), "horizontal diffusion (turbulent mixing) [on/off]")
    ("eqs.todo_mm.sedi", po::value<bool>(), "sedimentation [on/off]")
    //assumed initial aerosol parameters (for activation parametrisation only)
    ("eqs.todo_mm.mean_rd", po::value<string>(), "dry aerosol mean radii [m]") 
    ("eqs.todo_mm.sdev_rd", po::value<string>(), "geometric standard deviation [1]") 
    ("eqs.todo_mm.n_tot", po::value<string>(), "total concentration [m-3]") 
    ("eqs.todo_mm.chem_b", po::value<string>(), "chemical composition parameter [1]")

    ("eqs.todo_sdm.adve", po::value<bool>(), "advection [on/off]")
    ("eqs.todo_sdm.cond", po::value<bool>(), "condensation/evaporation [on/off]")
    ("eqs.todo_sdm.sedi", po::value<bool>(), "sedimentation [on/off]")
    ("eqs.todo_sdm.chem", po::value<bool>(), "chemistry [on/off]")
    ("eqs.todo_sdm.coal", po::value<bool>(), "coalescence [on/off]")
    ("eqs.todo_sdm.xi", po::value<string>()->default_value("p2"), "definition of xi (id, ln, p2, p3)")
    ("eqs.todo_sdm.adve.algo", po::value<string>()->default_value("euler"), "advection ODE solver type (euler, mmid, rk4)")
    ("eqs.todo_sdm.cond.algo", po::value<string>()->default_value("euler"), "condensation/evaporation ODE solver type (euler, mmid, rk4)")
    ("eqs.todo_sdm.sedi.algo", po::value<string>()->default_value("euler"), "sedimentation ODE solver type (euler, mmid, rk4)")
    ("eqs.todo_sdm.chem.algo", po::value<string>()->default_value("euler"), "chemistry ODE solver type (euler, mmid, rk4)")
    ("eqs.todo_sdm.adve.sstp", po::value<int>()->default_value(1), "number of substeps for advection")
    ("eqs.todo_sdm.cond.sstp", po::value<int>()->default_value(50), "number of substeps for condensation/evaporation")
    ("eqs.todo_sdm.sedi.sstp", po::value<int>()->default_value(1), "number of substeps for sedimentation")
    ("eqs.todo_sdm.chem.sstp", po::value<int>()->default_value(1), "number of substeps for chemistry")
    ("eqs.todo_sdm.coal.sstp", po::value<int>()->default_value(1), "number of substeps for coalescence")
    ("eqs.todo_sdm.sd_conc_mean", po::value<string>(), "mean super-droplet density per cell") 
    ("eqs.todo_sdm.min_rd", po::value<string>()->default_value("1e-9"), "minimum initial dry aerosol radius [m]") 
    ("eqs.todo_sdm.max_rd", po::value<string>()->default_value("1e-6"), "maximum initial dry aerosol radius [m]")
    // parameters for dry aerosol two-mode lognormal distribution
    ("eqs.todo_sdm.mean_rd1", po::value<string>(), "first mode dry aerosol mean radii [m]") 
    ("eqs.todo_sdm.mean_rd2", po::value<string>(), "second mode dry aerosol mean radii[m]")
    ("eqs.todo_sdm.sdev_rd1", po::value<string>(), "first mode geometric standard deviation [1]") 
    ("eqs.todo_sdm.sdev_rd2", po::value<string>(), "second mode geometric standard deviation [1]")
    ("eqs.todo_sdm.n1_tot", po::value<string>(), "first mode total concentration [m-3]") 
    ("eqs.todo_sdm.n2_tot", po::value<string>(), "second mode total concentration [m-3]")
    ("eqs.todo_sdm.kappa", po::value<string>(), "hygroscopicity parameter kappa [1]")
    ("eqs.todo_sdm.out_m0", po::value<string>()->default_value(".5e-6:25e-6"), "radius ranges for the 0-th moments [m]")
    ("eqs.todo_sdm.out_m1", po::value<string>()->default_value(".5e-6:25e-6"), "radius ranges for the 1-st moments [m]")
    ("eqs.todo_sdm.out_m2", po::value<string>()->default_value(".5e-6:25e-6"), "radius ranges for the 2-nd moments [m]")
    ("eqs.todo_sdm.out_m3", po::value<string>()->default_value(".5e-6:25e-6"), "radius ranges for the 3-rd moments [m]")
    ("eqs.todo_sdm.out_m6", po::value<string>()->default_value(""), "radius ranges for the 6-th moments [m]")
    // initial parameters of chemical compounds:
    // gas phase in the air
    ("eqs.todo_sdm.chem.env_SO2", po::value<string>()->default_value("1e-10"), "volume mixing ratio of SO2 [1]")
    ("eqs.todo_sdm.chem.env_O3", po::value<string>()->default_value("1e-10"), "volume mixing ratio of O3 [1]")
    ("eqs.todo_sdm.chem.env_H2O2", po::value<string>()->default_value("1e-10"), "volume mixing ratio of H2O2 [1]")
    // aqueous phase in the droplets
    //TODO initial conditions calculated to be consistent with environment and kappa!!!
    ("eqs.todo_sdm.chem.ini_c_H", po::value<string>()->default_value("1e-12"), "initial mass of H+ [kg]") //dissociation of pure water: 1e-4
    ("eqs.todo_sdm.chem.ini_c_OH", po::value<string>()->default_value("1e-25"), "initial mass of OH- [kg]") //dissociation of pure water:1e-4
    ("eqs.todo_sdm.chem.ini_c_SO2", po::value<string>()->default_value("0"), "initial mass of SO2*H2O [kg]")//1e-22
    ("eqs.todo_sdm.chem.ini_c_O3", po::value<string>()->default_value("0"), "initial mass of O3*H2O2 [kg]")
    ("eqs.todo_sdm.chem.ini_c_H2O2", po::value<string>()->default_value("0"), "initial mass of H2O2*H2O [kg]")
    ("eqs.todo_sdm.chem.ini_c_HSO3", po::value<string>()->default_value("1e-12"), "initial mass of HSO3- [kg]")
    ("eqs.todo_sdm.chem.ini_c_SO3", po::value<string>()->default_value("0"), "initial mass of SO3-- [kg]")
    ;
}

template <typename real_t>
eqs<real_t> *opt_eqs(
  const po::variables_map& vm, 
  const grd<real_t> &grid, 
  const ini<real_t> &intcond, 
  const vel<real_t> &velocity
) 
{
  string initype= vm.count("eqs") ? vm["eqs"].as<string>() : "<unspecified>";
  if (initype == "scalar_advection")
    return new eqs_scalar_advection<real_t>();
  else 
  if (initype == "isentropic")
  {
    if (!vm.count("eqs.isentropic.nlev")) error_macro("TODO")
    if (!vm.count("eqs.isentropic.abslev")) error_macro("TODO")
    return new eqs_isentropic<real_t>(grid, 
      vm["eqs.isentropic.nlev"].as<int>(),
      real_cast<real_t>(vm, "eqs.isentropic.p_top") * si::pascals,
      real_cast<real_t>(vm, "eqs.isentropic.theta_frst") * si::kelvins,
      vm["eqs.isentropic.abslev"].as<int>(),
      real_cast<real_t>(vm, "eqs.isentropic.absamp")
    );
  }
  else
  if (initype == "todo_bulk")
    return new eqs_todo_bulk_ode<real_t>(grid,
      map<enum eqs_todo_bulk<real_t>::processes, bool>({
        {eqs_todo_bulk<real_t>::cevp, vm["eqs.todo_bulk.cevp"].as<bool>()},
        {eqs_todo_bulk<real_t>::conv, vm["eqs.todo_bulk.conv"].as<bool>()},
        {eqs_todo_bulk<real_t>::clct, vm["eqs.todo_bulk.clct"].as<bool>()},
        {eqs_todo_bulk<real_t>::sedi, vm["eqs.todo_bulk.sedi"].as<bool>()},
        {eqs_todo_bulk<real_t>::revp, vm["eqs.todo_bulk.revp"].as<bool>()}
      })
    );
  else
  if (initype == "todo_mm")
    return new eqs_todo_mm<real_t>(grid, 
       map<enum eqs_todo_mm<real_t>::processes, bool>({
        {eqs_todo_mm<real_t>::act,   vm["eqs.todo_mm.act"].as<bool>()},
        {eqs_todo_mm<real_t>::cond,  vm["eqs.todo_mm.cond"].as<bool>()},
        {eqs_todo_mm<real_t>::acc,   vm["eqs.todo_mm.acc"].as<bool>()},
        {eqs_todo_mm<real_t>::autoc, vm["eqs.todo_mm.autoc"].as<bool>()},
        {eqs_todo_mm<real_t>::turb,  vm["eqs.todo_mm.turb"].as<bool>()},
        {eqs_todo_mm<real_t>::sedi,  vm["eqs.todo_mm.sedi"].as<bool>()}
      }),
      real_cast<real_t>(vm, "eqs.todo_mm.mean_rd"), 
      real_cast<real_t>(vm, "eqs.todo_mm.sdev_rd"), 
      real_cast<real_t>(vm, "eqs.todo_mm.n_tot"), 
      real_cast<real_t>(vm, "eqs.todo_mm.chem_b")
    );
  else 
  if (initype == "todo_sdm")
  {
    map<string, enum sdm::ode_algos> map_algo({
      {"euler", sdm::euler},
      {"mmid",  sdm::mmid},
      {"rk4",   sdm::rk4}
    });
    map<string, enum sdm::xi_dfntns> map_xid({
      {"id", sdm::id},
      {"ln", sdm::ln},
      {"p2", sdm::p2},
      {"p3", sdm::p3}
    });

    // parsing moment list
    map<int, vector<pair<
      quantity<si::length, real_t>, 
      quantity<si::length, real_t>
    >>> outmoments;
    for (int &k : std::list<int>({0,1,2,3,6})) 
    {
      ostringstream opt;
      opt << "eqs.todo_sdm.out_m" << k;
      string str = vm[opt.str()].as<string>();

      std::vector<pair<std::string, std::string>> ranges;
      std::string::iterator first = str.begin();
      std::string::iterator last  = str.end();

      const bool result = boost::spirit::qi::phrase_parse(first, last, 
        *( 
          *(boost::spirit::qi::char_-":")  >> 
          boost::spirit::qi::lit(":") >> 
          *(boost::spirit::qi::char_-";") >> 
          -boost::spirit::qi::lit(";") 
        ),
        boost::spirit::ascii::space, ranges
      );    
      if (!result || first != last) error_macro("failed to parse string: " << str)
      for (auto &range : ranges)
      {
        try 
        {
          outmoments[k].push_back({
            boost::lexical_cast<real_t>(range.first) * si::metres,
            boost::lexical_cast<real_t>(range.second) * si::metres
          });
        }
        catch (boost::bad_lexical_cast &)  
        {
          error_macro("failed to convert " << range.first << " ... " << range.second
            << " into floating point value range"
          )
        }
      }
    }

    // calling the ctor
    return new eqs_todo_sdm<real_t>(grid, velocity,
      map<enum sdm::processes, bool>({
        {sdm::adve, vm["eqs.todo_sdm.adve"].as<bool>()},
        {sdm::cond, vm["eqs.todo_sdm.cond"].as<bool>()},
        {sdm::sedi, vm["eqs.todo_sdm.sedi"].as<bool>()},
        {sdm::coal, vm["eqs.todo_sdm.coal"].as<bool>()},
        {sdm::chem, vm["eqs.todo_sdm.chem"].as<bool>()},
      }),
      map_xid.at(vm["eqs.todo_sdm.xi"].as<string>()),
      map_algo.at(vm["eqs.todo_sdm.adve.algo"].as<string>()),
      map_algo.at(vm["eqs.todo_sdm.sedi.algo"].as<string>()),
      map_algo.at(vm["eqs.todo_sdm.cond.algo"].as<string>()),
      map_algo.at(vm["eqs.todo_sdm.chem.algo"].as<string>()),
      vm["eqs.todo_sdm.adve.sstp"].as<int>(),
      vm["eqs.todo_sdm.sedi.sstp"].as<int>(),
      vm["eqs.todo_sdm.cond.sstp"].as<int>(),
      vm["eqs.todo_sdm.chem.sstp"].as<int>(),
      vm["eqs.todo_sdm.coal.sstp"].as<int>(),
      real_cast<real_t>(vm, "eqs.todo_sdm.sd_conc_mean"),
      real_cast<real_t>(vm, "eqs.todo_sdm.min_rd"),
      real_cast<real_t>(vm, "eqs.todo_sdm.max_rd"),
      real_cast<real_t>(vm, "eqs.todo_sdm.mean_rd1"), 
      real_cast<real_t>(vm, "eqs.todo_sdm.mean_rd2"),
      real_cast<real_t>(vm, "eqs.todo_sdm.sdev_rd1"), 
      real_cast<real_t>(vm, "eqs.todo_sdm.sdev_rd2"),
      real_cast<real_t>(vm, "eqs.todo_sdm.n1_tot"), 
      real_cast<real_t>(vm, "eqs.todo_sdm.n2_tot"),
      real_cast<real_t>(vm, "eqs.todo_sdm.kappa"),
      map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> ({
        {sdm::gSO2,  real_cast<real_t>(vm, "eqs.todo_sdm.chem.env_SO2")},
        {sdm::gO3,   real_cast<real_t>(vm, "eqs.todo_sdm.chem.env_O3")},
        {sdm::gH2O2, real_cast<real_t>(vm, "eqs.todo_sdm.chem.env_H2O2")},
      }),
      map<enum sdm::chem_aq, quantity<si::mass, real_t>> ({
        {sdm::H,    real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_H")   *si::kilograms},
        {sdm::OH,    real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_OH") *si::kilograms},
        {sdm::SO2,  real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_SO2") *si::kilograms},
        {sdm::O3,   real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_O3")  *si::kilograms},
        {sdm::H2O2, real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_H2O2")*si::kilograms},
        {sdm::HSO3, real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_HSO3")*si::kilograms},
        {sdm::SO3,  real_cast<real_t>(vm, "eqs.todo_sdm.chem.ini_c_SO3") *si::kilograms}
      }),
      outmoments
    );
  }
  else 
  error_macro("unsupported equation system: " << initype)
}
