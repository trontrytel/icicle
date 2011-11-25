/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ARR_HPP
#  define ARR_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting, ...

template <class unit, typename real_t>
class arr : root
{
  private: Array<quantity<unit, real_t>, 3> *arr_ijk, *arr_jki, *arr_kij;

  public: arr(Range rng_x, Range rng_y, Range rng_z)
  {
    arr_ijk = new Array<quantity<unit, real_t>, 3>(rng_x, rng_y, rng_z);
    arr_jki = array_view(arr_ijk, secondRank, thirdRank, firstRank);
    arr_kij = array_view(arr_ijk, thirdRank, firstRank, secondRank);
  }

  public: ~arr()
  {
    delete arr_kij;
    delete arr_jki;
    delete arr_ijk;
  }

  private: Array<quantity<unit, real_t>, 3>*array_view(
    Array<quantity<unit, real_t>, 3>* a, int d1, int d2, int d3
  )
  {
    GeneralArrayStorage<3> storage;
    storage.ordering() = d1, d2, d3; 
    storage.base() = a->lbound(d1), a->lbound(d2), a->lbound(d3);
    return new Array<quantity<unit, real_t>, 3>( 
      a->dataFirst(), 
      TinyVector<int, 3>(a->extent(d1), a->extent(d2), a->extent(d3)), 
      storage
    );  
  }

  public: Array<quantity<unit, real_t>, 3> &ijk() { return *arr_ijk; }
  public: Array<quantity<unit, real_t>, 3> &jki() { return *arr_jki; }
  public: Array<quantity<unit, real_t>, 3> &kij() { return *arr_kij; }
};
#endif
