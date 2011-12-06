  public: void op_ijk(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz 
  )   
  {
    op<idx_ijk>(0, psi, tmp, i, j, k, n, step, Cx, Cy, Cz);
  }

  public: void op_jki(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz 
  )   
  {
    op<idx_jki>(1, psi, tmp, j, k, i, n, step, Cx, Cy, Cz);
  }

  public: void op_kij(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz 
  )   
  {
    op<idx_kij>(2, psi, tmp, k, i, j, n, step, Cx, Cy, Cz);
  }
