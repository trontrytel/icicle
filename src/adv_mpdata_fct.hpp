/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_mpdata_fct class with an implementation of the flux-corrected MPDATA scheme
 */
#ifndef ADV_MPDATA_FCT_HPP
#  define ADV_MPDATA_FCT_HPP

#  include "adv_mpdata.hpp"

/** @brief
 *  Flux Corrected Transport (aka non-oscillatory, monotonic, sign-preserving)
 *  option for MPDATA (for the transport of positive scalars only)
 */
template <typename real_t> 
class adv_mpdata_fct : public adv_mpdata<real_t> 
{
  public: const int stencil_extent() 
  {
    int halo = (adv_mpdata<real_t>::stencil_extent() - 1) / 2; 
    halo += 1; // cf. ii, jj & kk below
    return 1 + 2 * halo;
  }
  public: const int num_vctr_caches() 
  { 
    return 1 + (iord < 3 ? 3 : 6); 
  }
  public: const int num_sclr_caches() { return 2; } 
#  define psi_min (*tmp_s[0])
#  define psi_max (*tmp_s[1])

  private: int iord;
  private: bool cross_terms, third_order;
  private: grd_arakawa_c_lorenz<real_t> *grid;
  public: adv_mpdata_fct(grd_arakawa_c_lorenz<real_t> *grid, int iord, bool cross_terms, bool third_order) : 
    adv_mpdata<real_t>(grid, iord, cross_terms, third_order), 
    iord(iord), cross_terms(cross_terms), third_order(third_order), grid(grid)
  {
    warning_macro("enabling FCT corrections reduces the order of the scheme")
  }

  public: 
  template <class mpdata_op3D_class>
  class op3D : public adv<real_t>::op3D
  {
    private:
    template <class idx>
    class indices 
    { 
      public: mtx::rng iv, ii, jj, kk;
      public: idx ii_jj_kk, iim1_jj_kk, iip1_jj_kk, ii_jjm1_kk, ii_jjp1_kk, ii_jj_kkm1, ii_jj_kkp1;
      public: idx iv_j_k, ivph_j_k, ivmh_j_k, iv_jph_k, iv_jmh_k, iv_j_kph, iv_j_kmh,
        ivp1mh_j_k, ivp1ph_j_k, ivp2_j_k, ivp1_jmh_k, ivp1_jm1_k, ivp1_jph_k, ivp1_jp1_k,
        ivp1_j_kmh, ivp1_j_km1, ivp1_j_kph, ivp1_j_kp1, ivp1_j_k, ivm1_j_k, iv_jm1_k, iv_jp1_k, 
        iv_j_km1, iv_j_kp1;
      public: indices(
        const mtx::rng &i, 
        const mtx::rng &j, 
        const mtx::rng &k, 
        const grd_arakawa_c_lorenz<real_t> &grid,
        int cross_terms,
        int iord
      ) : 
        ii(i.first() - 1, i.last() + 1),
        jj(j.first() - 1, j.last() + 1),
        kk(k.first() - 1, k.last() + 1),

        ii_jj_kk(ii, jj, kk),
        iim1_jj_kk(ii-1, jj, kk), 
        iip1_jj_kk(ii+1, jj, kk), 
        ii_jjm1_kk(ii, jj-1, kk), 
        ii_jjp1_kk(ii, jj+1, kk), 
        ii_jj_kkm1(ii, jj, kk-1), 
        ii_jj_kkp1(ii, jj, kk+1),

        iv(i.first()-1, i.last()), // as in mpdata_U, we compute u_{i+1/2} for iv=(i-1, ... i) instead of u_{i+1/2} and u_{i-1/2} for all i

        iv_j_k(iv, j, k), 
        ivph_j_k(iv + grid.p_half, j, k), 
        ivmh_j_k(iv - grid.m_half, j, k), 
        iv_jph_k(iv, j + grid.p_half, k),
        iv_jmh_k(iv, j - grid.m_half, k), 
        iv_j_kph(iv, j, k + grid.p_half), 
        iv_j_kmh(iv, j, k - grid.m_half),
        ivp1mh_j_k(iv + 1 - grid.m_half, j, k),
        ivp1ph_j_k(iv + 1 + grid.p_half, j, k),
        ivp2_j_k(iv + 2, j, k),
        ivp1_jmh_k(iv + 1, j - grid.m_half, k),
        ivp1_jm1_k(iv + 1, j - 1, k),
        ivp1_jph_k(iv + 1, j + grid.p_half, k),
        ivp1_jp1_k(iv + 1, j + 1, k),
        ivp1_j_kmh(iv + 1, j, k - grid.m_half),
        ivp1_j_km1(iv + 1, j, k - 1),
        ivp1_j_kph(iv + 1, j, k + grid.p_half),
        ivp1_j_kp1(iv + 1, j, k + 1),
        ivp1_j_k(iv + 1, j, k),
        ivm1_j_k(iv - 1, j, k), 
        iv_jm1_k(iv, j - 1, k), 
        iv_jp1_k(iv, j + 1, k), 
        iv_j_km1(iv, j, k - 1), 
        iv_j_kp1(iv, j, k + 1)
      { }
    };  

    // helper adv ops
    private: mpdata_op3D_class mpdata;
    private: typename adv_upstream<real_t>::op3D upstream;

    // private members 
    private: indices<mtx::idx_ijk> idxx;
    private: indices<mtx::idx_jki> idxy;
    private: indices<mtx::idx_kij> idxz;
    private: mtx::arr<real_t> **tmp_s, **tmp_v;
    private: int iord;
    private: bool positive_definite;

    // ctor
    public: op3D(
      const mtx::idx &ijk,
      const grd_arakawa_c_lorenz<real_t> &grid,
      mtx::arr<real_t> **tmp_s,
      mtx::arr<real_t> **tmp_v,
      int cross_terms, int iord, int third_order, bool positive_definite
    ) :
      adv<real_t>::op3D(ijk),
      mpdata(mtx::idx_ijk(
        mtx::rng(ijk.i.first() - 1, ijk.i.last() + 1),
        mtx::rng(ijk.j.first() - 1, ijk.j.last() + 1),
        mtx::rng(ijk.k.first() - 1, ijk.k.last() + 1)
      ), grid, tmp_s, tmp_v, cross_terms, iord, third_order),
      upstream(ijk, grid),
      tmp_s(tmp_s), tmp_v(tmp_v), iord(iord), positive_definite(positive_definite),
      idxx(ijk.i, ijk.j, ijk.k, grid, cross_terms, iord),
      idxy(ijk.j, ijk.k, ijk.i, grid, cross_terms, iord),
      idxz(ijk.k, ijk.i, ijk.j, grid, cross_terms, iord)
    { }

    public: virtual void operator()(
      mtx::arr<real_t> *psi[],
      const int n,
      const int step,
      const mtx::arr<real_t> * const Cx,
      const mtx::arr<real_t> * const Cy,
      const mtx::arr<real_t> * const Cz
    )
    {

      assert(isfinite(sum((*psi[n])(idxx.ii_jj_kk))));
      assert(isfinite(sum((*psi[n])(idxx.iim1_jj_kk))));
      assert(isfinite(sum((*psi[n])(idxx.iip1_jj_kk))));
      assert(isfinite(sum((*psi[n])(idxx.ii_jjm1_kk))));
      assert(isfinite(sum((*psi[n])(idxx.ii_jjp1_kk))));
      assert(isfinite(sum((*psi[n])(idxx.ii_jj_kkm1))));
      assert(isfinite(sum((*psi[n])(idxx.ii_jj_kkp1))));

#  define mpdata_fct_minmax(fun, psi_, n_, idx) fun( \
     (*psi_[n_])(idx.ii_jj_kk), fun( \
     (*psi_[n_])(idx.iim1_jj_kk), fun( \
     (*psi_[n_])(idx.iip1_jj_kk), fun( \
     (*psi_[n_])(idx.ii_jjm1_kk), fun( \
     (*psi_[n_])(idx.ii_jjp1_kk), fun( \
     (*psi_[n_])(idx.ii_jj_kkm1), \
     (*psi_[n_])(idx.ii_jj_kkp1)))))) \
   ) 

    /// \f$ \psi^{max}_{i}=max_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n
    /// \f$ \psi^{min}_{i}=min_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n
    /// eq.(20a, 20b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
      if (step == 1)
      {
        // calculating psi_min and psi_max from the previous time step
        psi_min(idxx.ii_jj_kk) = mpdata_fct_minmax(min, psi, n, idxx);
        psi_max(idxx.ii_jj_kk) = mpdata_fct_minmax(max, psi, n, idxx);
        // performing standard upstream advection
        upstream(psi, n, 1, Cx, Cy, Cz);
      }
      else
      {
        // calculating psi_min and psi_max from the previous time step and previous iord
        psi_min(idxx.ii_jj_kk) = min(psi_min(idxx.ii_jj_kk), mpdata_fct_minmax(min, psi, n, idxx));
        psi_max(idxx.ii_jj_kk) = max(psi_max(idxx.ii_jj_kk), mpdata_fct_minmax(max, psi, n, idxx));

        // calculating Cx_mon, Cy_mon, Cz_mon
        int 
          x_old = -1, y_old = -1, z_old = -1, 
          x_new =  1, y_new =  2, z_new =  3; 
        if (iord > 2)
        {   
          if (step > 2)
          {
            x_old = 1, y_old = 2, z_old = 3;
            if (step % 2 == 0)  
              x_old += 3, y_old += 3, z_old += 3;
            else 
              x_new += 3, y_new += 3, z_new += 3;
          }
        }
  
        const mtx::arr<real_t> 
          * const Cx_unco = (step < 3 ? Cx : tmp_v[x_old]), 
          *       Cx_corr = tmp_v[x_new],
          *       Cx_mono = tmp_v[0],
          * const Cy_unco = (step < 3 ? Cy : tmp_v[y_old]), 
          *       Cy_corr = tmp_v[y_new],
          *       Cy_mono = tmp_v[0],
          * const Cz_unco = (step < 3 ? Cz : tmp_v[z_old]), 
          *       Cz_corr = tmp_v[z_new],
          *       Cz_mono = tmp_v[0];
  
        if (this->do_x())
          mpdata.mpdata_U(Cx_corr, psi, n, step, mpdata.indcs_x, *Cx_unco, *Cy_unco, *Cz_unco);
        if (this->do_y())
          mpdata.mpdata_U(Cy_corr, psi, n, step, mpdata.indcs_y, *Cy_unco, *Cz_unco, *Cx_unco);
        if (this->do_z())
          mpdata.mpdata_U(Cz_corr, psi, n, step, mpdata.indcs_z, *Cz_unco, *Cx_unco, *Cy_unco); 
   
        // performing upstream advection using the ''monotonic'' velocities (logic from adv::op3D)
        *psi[n+1] = *psi[0]; // TODO: at least this should be placed in adv... and the leapfrog & upstream in adv_dimsplit?
        if (this->do_x())
        {
          fct_helper(psi, *Cx_mono, *Cx_corr, *Cy_corr, *Cz_corr, idxx, n);
          upstream.op1D(psi, upstream.indcs_x, n, 1, Cx_mono, NULL, NULL); 
        }
        if (this->do_y())
        {
          fct_helper(psi, *Cy_mono, *Cy_corr, *Cz_corr, *Cx_corr, idxy, n); 
          upstream.op1D(psi, upstream.indcs_y, n, 1, Cy_mono, NULL, NULL); 
        }
        if (this->do_z())
        {
          fct_helper(psi, *Cz_mono, *Cz_corr, *Cx_corr, *Cy_corr, idxz, n); 
          upstream.op1D(psi, upstream.indcs_z, n, 1, Cz_mono, NULL, NULL); 
        }
      }
    }

    private:
    template <class indices>
    void fct_helper(
      mtx::arr<real_t> *psi[], 
      const mtx::arr<real_t> &C_mon_x, 
      const mtx::arr<real_t> &C_adf_x, const mtx::arr<real_t> &C_adf_y, const mtx::arr<real_t> &C_adf_z, 
      const indices &idx,
      int n
    )  
    {
///
/// \f$ \beta^{\uparrow}_{i} = \frac { \psi^{max}_{i}- \psi^{*}_{i} }
/// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i-1/2}]^{+} \psi^{*}_{i-1} - 
/// [u^{I}_{i+1/2}]^{-} \psi^{*}_{i+1} \right)  } \f$ \n
/// \f$ \beta^{\downarrow}_{i} = \frac { \psi^{*}_{i}- \psi^{min}_{i} }
/// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i+1/2}]^{+} \psi^{*}_{i} - 
/// [u^{I}_{i-1/2}]^{-} \psi^{*}_{i} \right)  } \f$ \n
/// eq.(19a, 19b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

#  define mpdata_fct_beta_up(i_j_k, imh_j_k, iph_j_k, i_jmh_k, i_jph_k, i_j_kmh, i_j_kph, im1_j_k, ip1_j_k, i_jm1_k, i_jp1_k, i_j_km1, i_j_kp1) adv_mpdata<real_t>::frac( \
     psi_max(i_j_k) - (*psi[n])(i_j_k), \
     adv<real_t>::pospart(C_adf_x(imh_j_k)) * (*psi[n])(im1_j_k) - \
     adv<real_t>::negpart(C_adf_x(iph_j_k)) * (*psi[n])(ip1_j_k) + \
     adv<real_t>::pospart(C_adf_y(i_jmh_k)) * (*psi[n])(i_jm1_k) - \
     adv<real_t>::negpart(C_adf_y(i_jph_k)) * (*psi[n])(i_jp1_k) + \
     adv<real_t>::pospart(C_adf_z(i_j_kmh)) * (*psi[n])(i_j_km1) - \
     adv<real_t>::negpart(C_adf_z(i_j_kph)) * (*psi[n])(i_j_kp1)   \
   )
#  define mpdata_fct_beta_dn(i_j_k, iph_j_k, imh_j_k, i_jph_k, i_jmh_k, i_j_kph, i_j_kmh) adv_mpdata<real_t>::frac(\
     (*psi[n])(i_j_k) - psi_min(i_j_k), \
     adv<real_t>::pospart(C_adf_x(iph_j_k)) * (*psi[n])(i_j_k) - \
     adv<real_t>::negpart(C_adf_x(imh_j_k)) * (*psi[n])(i_j_k) + \
     adv<real_t>::pospart(C_adf_y(i_jph_k)) * (*psi[n])(i_j_k) - \
     adv<real_t>::negpart(C_adf_y(i_jmh_k)) * (*psi[n])(i_j_k) + \
     adv<real_t>::pospart(C_adf_z(i_j_kph)) * (*psi[n])(i_j_k) - \
     adv<real_t>::negpart(C_adf_z(i_j_kmh)) * (*psi[n])(i_j_k)   \
   )
      /// nonoscillatory antidiffusive velocity: \n
      /// \f$ U^{MON}_{i+1/2}=min(1,\beta ^{\downarrow}_i,\beta ^{\uparrow} _{i+1})[U_{i+1/2}]^{+} 
      /// + min(1,\beta^{\uparrow}_{i},\beta^{\downarrow}_{i+1/2})[u_{i+1/2}]^{-} \f$ \n
      /// where \f$ [\cdot]^{+}=max(\cdot,0) \f$ and \f$ [\cdot]^{-}=min(\cdot,0) \f$ \n
      /// eq.(18) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

/*
cerr << "aqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.iv_j_k))));
cerr << "bqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.ivp1_j_k))));
cerr << "cqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.ivm1_j_k))));
cerr << "dqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.ivp2_j_k))));
cerr << "eqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.iv_jm1_k))));
cerr << "fqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.ivp1_jp1_k))));
cerr << "gqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.iv_j_km1))));
cerr << "hqq" << endl;
    assert(isfinite(sum((*psi[n])(idx.ivp1_j_kp1))));

cerr << "iqq" << endl;
cerr << C_adf_x(idx.ivph_j_k) << endl;
    assert(isfinite(sum(C_adf_x(idx.ivph_j_k))));
cerr << "jqq" << endl;
    assert(isfinite(sum(C_adf_x(idx.ivmh_j_k))));
cerr << "kqq" << endl;
    assert(isfinite(sum(C_adf_x(idx.ivp1ph_j_k))));
cerr << "lqq" << endl;
    assert(isfinite(sum(C_adf_x(idx.ivp1mh_j_k))));

cerr << "mqq" << endl;
    assert(isfinite(sum(C_adf_y(idx.iv_jph_k))));
cerr << "nqq" << endl;
    assert(isfinite(sum(C_adf_y(idx.iv_jmh_k))));
cerr << "oqq" << endl;
    assert(isfinite(sum(C_adf_y(idx.ivp1_jph_k))));
cerr << "pqq" << endl;
    assert(isfinite(sum(C_adf_y(idx.ivp1_jmh_k))));

cerr << "qqq" << endl;
    assert(isfinite(sum(C_adf_z(idx.iv_j_kph))));
cerr << "rqq" << endl;
    assert(isfinite(sum(C_adf_z(idx.iv_j_kmh))));
cerr << "sqq" << endl;
    assert(isfinite(sum(C_adf_z(idx.ivp1_j_kph))));
cerr << "tqq" << endl;
    assert(isfinite(sum(C_adf_z(idx.ivp1_j_kmh))));
*/

      if (positive_definite) 
        C_mon_x(idx.ivph_j_k) = C_adf_x(idx.ivph_j_k) * where(
          C_adf_x(idx.ivph_j_k) > 0,
          min(1, min(
            mpdata_fct_beta_dn(idx.iv_j_k, idx.ivph_j_k, idx.ivmh_j_k, idx.iv_jph_k, idx.iv_jmh_k, idx.iv_j_kph, idx.iv_j_kmh), 
            mpdata_fct_beta_up(idx.ivp1_j_k, idx.ivp1mh_j_k, idx.ivp1ph_j_k, idx.ivp1_jmh_k, idx.ivp1_jph_k, idx.ivp1_j_kmh, idx.ivp1_j_kph, idx.iv_j_k, idx.ivp2_j_k, idx.ivp1_jm1_k, idx.ivp1_jp1_k, idx.ivp1_j_km1, idx.ivp1_j_kp1))
          ),
          min(1, min(
            mpdata_fct_beta_up(idx.iv_j_k, idx.ivmh_j_k, idx.ivph_j_k, idx.iv_jmh_k, idx.iv_jph_k, idx.iv_j_kmh, idx.iv_j_kph, idx.ivm1_j_k, idx.ivp1_j_k, idx.iv_jm1_k, idx.iv_jp1_k, idx.iv_j_km1, idx.iv_j_kp1), 
            mpdata_fct_beta_dn(idx.ivp1_j_k, idx.ivp1ph_j_k, idx.ivp1mh_j_k, idx.ivp1_jph_k, idx.ivp1_jmh_k, idx.ivp1_j_kph, idx.ivp1_j_kmh))
          )
        );
      else
        C_mon_x(idx.ivph_j_k) = C_adf_x(idx.ivph_j_k) * 
          where(
            C_adf_x(idx.ivph_j_k) > 0,
            where(
              (*psi[n])(idx.iv_j_k) > 0,
              min(1, min(
                mpdata_fct_beta_dn(idx.iv_j_k, idx.ivph_j_k, idx.ivmh_j_k, idx.iv_jph_k, idx.iv_jmh_k, idx.iv_j_kph, idx.iv_j_kmh), 
                mpdata_fct_beta_up(idx.ivp1_j_k, idx.ivp1mh_j_k, idx.ivp1ph_j_k, idx.ivp1_jmh_k, idx.ivp1_jph_k, idx.ivp1_j_kmh, idx.ivp1_j_kph, idx.iv_j_k, idx.ivp2_j_k, idx.ivp1_jm1_k, idx.ivp1_jp1_k, idx.ivp1_j_km1, idx.ivp1_j_kp1)
              )),
              min(1, min(
                mpdata_fct_beta_up(idx.iv_j_k, idx.ivmh_j_k, idx.ivph_j_k, idx.iv_jmh_k, idx.iv_jph_k, idx.iv_j_kmh, idx.iv_j_kph, idx.ivm1_j_k, idx.ivp1_j_k, idx.iv_jm1_k, idx.iv_jp1_k, idx.iv_j_km1, idx.iv_j_kp1), 
                mpdata_fct_beta_dn(idx.ivp1_j_k, idx.ivp1ph_j_k, idx.ivp1mh_j_k, idx.ivp1_jph_k, idx.ivp1_jmh_k, idx.ivp1_j_kph, idx.ivp1_j_kmh)
              ))
            ),
            where(
              (*psi[n])(idx.ivp1_j_k) > 0,
              min(1, min(
                mpdata_fct_beta_up(idx.iv_j_k, idx.ivmh_j_k, idx.ivph_j_k, idx.iv_jmh_k, idx.iv_jph_k, idx.iv_j_kmh, idx.iv_j_kph, idx.ivm1_j_k, idx.ivp1_j_k, idx.iv_jm1_k, idx.iv_jp1_k, idx.iv_j_km1, idx.iv_j_kp1), 
                mpdata_fct_beta_dn(idx.ivp1_j_k, idx.ivp1ph_j_k, idx.ivp1mh_j_k, idx.ivp1_jph_k, idx.ivp1_jmh_k, idx.ivp1_j_kph, idx.ivp1_j_kmh)
              )),
              min(1, min(
                mpdata_fct_beta_dn(idx.iv_j_k, idx.ivph_j_k, idx.ivmh_j_k, idx.iv_jph_k, idx.iv_jmh_k, idx.iv_j_kph, idx.iv_j_kmh), 
                mpdata_fct_beta_up(idx.ivp1_j_k, idx.ivp1mh_j_k, idx.ivp1ph_j_k, idx.ivp1_jmh_k, idx.ivp1_jph_k, idx.ivp1_j_kmh, idx.ivp1_j_kph, idx.iv_j_k, idx.ivp2_j_k, idx.ivp1_jm1_k, idx.ivp1_jp1_k, idx.ivp1_j_km1, idx.ivp1_j_kp1)
              ))
            )
          );
    }
  };
  
  public: typename adv<real_t>::op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **tmp_s,
    mtx::arr<real_t> **tmp_v,
    bool positive_definite
  )
  {
    if (positive_definite)
      return new op3D<typename adv_mpdata<real_t>::template op3D<typename adv_mpdata<real_t>::aon_nil>>(
        ijk, *grid, tmp_s, tmp_v, cross_terms, iord, third_order, positive_definite
      );
    else
      return new op3D<typename adv_mpdata<real_t>::template op3D<typename adv_mpdata<real_t>::aon_abs>>(
        ijk, *grid, tmp_s, tmp_v, cross_terms, iord, third_order, positive_definite
      );
  }
};
#endif
