  public: void op_ijk(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz 
  )   
  {
    op<idx_ijk>(0, psi, tmp_s, tmp_v, i, j, k, n, step, Cx, Cy, Cz);
  }

  public: void op_jki(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz 
  )   
  {
    op<idx_jki>(1, psi, tmp_s, tmp_v, j, k, i, n, step, Cy, Cz, Cx);
  }

  public: void op_kij(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz 
  )   
  {
    op<idx_kij>(2, psi, tmp_s, tmp_v, k, i, j, n, step, Cz, Cx, Cy);
  }
