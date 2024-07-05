#include "GMRES.hpp"

void GMRES::solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
                  double epsilon) const {
  double d = d0;
  RMatrix H(M + 1, M);
  RVector e(M + 1), y(M), cs(M), sn(M);
  LazyVectors v(M + 1, r);
  Vector c(r), t(r);
  int iter = 0;
  for (; iter < max_iter;) {
    if ((iter >= min_iter) && (d < epsilon)) break;
    vout(1) << iteration;
    v[0] = (1 / d) * r;
    e[0] = d;
    for (int i = 1; i < M; ++i)
      e[i] = 0;
    int k = 0;
    for (; (k < M) && (iter < max_iter); ++k) {
      vout(10) << Date() << endl;
      c = B * v[k];
      t = A * c;
      vout(2) << k << " c " << c.ConsistentNorm() << " Ac " << norm(t) << endl;
      v[k + 1] = t;
      for (int i = 0; i <= k; ++i) {
        H[i][k] = t * v[i];
        v[k + 1] -= H[i][k] * v[i];
      }
      H[k + 1][k] = norm(v[k + 1]);
      vout(3) << k << " |v_k+1| " << H[k + 1][k] << endl;
      v[k + 1] *= (1.0 / H[k + 1][k]);
      for (int i = 0; i < k; ++i) {
        Scalar giv_r = H[i][k];
        Scalar giv_h = H[i + 1][k];
        H[i][k] = sn[i] * giv_r + cs[i] * giv_h;
        H[i + 1][k] = cs[i] * giv_r - sn[i] * giv_h;
      }
      Scalar giv_r = H[k][k];
      Scalar giv_h = H[k + 1][k];
      Scalar co, si;
      d = sqrt(giv_r * giv_r + giv_h * giv_h);
      vout(3) << k << " d " << d << endl;
      if (giv_r > 0) {
        co = giv_h / d;
        si = giv_r / d;
        H[k][k] = d;
      } else {
        co = -giv_h / d;
        si = -giv_r / d;
        H[k][k] = -d;
      }
      H[k + 1][k] = 0;
      giv_r = e[k];
      e[k] = si * giv_r;
      e[k + 1] = co * giv_r;
      cs[k] = co;
      sn[k] = si;
      ++iter;
      iteration.push_back(abs(e[k + 1]));
      if ((iter >= min_iter) && (abs(e[k + 1]) < epsilon)) {
        k++;
        break;
      }
      if (iter % printSteps == 0) { vout(1) << iteration; }
    }
    for (int i = k - 1; i >= 0; --i) {
      Scalar s = e[i];
      for (int j = i + 1; j < k; ++j)
        s -= H[i][j] * y[j];
      y[i] = s / H[i][i];
    }
    t.Clear();
    for (int i = 0; i < k; ++i)
      t += y[i] * v[i];
    c = B * t;
    u += c;
    r -= A * c;
    d = norm(r);
    iteration.replace_back(d);
  }
}

void GMRES::solve(const Operator &A, const Operator &B, Vectors &u, Vectors &r, double d0,
                  double eps) const {
  mout.StartBlock("GMRES");
  int size = int(u.size()); // solve "size" rhs
  vector<int> N(size);
  for (int n = 0; n < size; ++n)
    N[n] = n;
  Vectors U(u);
  Vectors R(r);
  const int M = 100;
  vector<double> dvec = R.norm();
  double d_0 = *max_element(dvec.begin(), dvec.begin() + size);
  double d = d_0;
  vector<double> epsilon(size);
  for (int n = 0; n < size; ++n)
    epsilon[n] = max(defaultEpsilon, defaultReduction * dvec[n]);
  double e = max(defaultEpsilon, defaultReduction * d);
  int iter_min = -1;
  double d_min = 1e20;
  double d_max = -1;

  Vectors c(U), t(U);
  int iter = 0;
  for (; iter < max_iter;) {
    for (int n = 0; n < size; ++n)
      if (dvec[n] < epsilon[n]) { // one rhs is completed; decreasing size
        if (iter_min == -1) iter_min = iter;
        u[N[n]] = U[n];
        r[N[n]] = R[n];
        U.erase(n);
        N.erase(N.begin() + n);
        d_min = min(d_min, dvec[n]);
        d_max = max(d_max, dvec[n]);
        dvec.erase(dvec.begin() + n);
        epsilon.erase(epsilon.begin() + n);
        c.erase(n);
        t.erase(n);
        --size;
        --n;
      }
    if (size == 0) break;
    vout(1) << "d(" << iter << ")= " << d << " (maximum of " << size << " rhs)" << endl;
    vector<Vectors> v(1, U);
    vector<vector<vector<Scalar>>> H(size);
    for (int n = 0; n < size; ++n) {
      H[n].resize(M + 1);
      for (int m = 0; m < M + 1; ++m)
        H[n][m].resize(M);
    }
    vector<vector<Scalar>> e(size), y(size), cs(size), sn(size);
    for (int n = 0; n < size; ++n) {
      e[n].resize(M + 1);
      y[n].resize(M);
      cs[n].resize(M);
      sn[n].resize(M);
    }
    for (int n = 0; n < size; ++n) {
      v[0][n] = 1 / dvec[n] * R[n];
      e[n][0] = dvec[n];
      for (int i = 1; i < M; ++i)
        e[n][i] = 0;
    }
    int k = 0;
    Vectors C(c), T(t);
    int ksize = size;
    vector<int> K(size);
    for (int n = 0; n < size; ++n)
      K[n] = n;
    for (; (k < M) && (iter < max_iter); ++k) {
      C = B * v[k];
      T = A * C;
      v.push_back(T);
      //                     v[k+1] = T;
      for (int n = 0; n < ksize; ++n) {
        for (int i = 0; i <= k; ++i) {
          H[n][i][k] = T[n] * v[i][n];
          v[k + 1][n] -= H[n][i][k] * v[i][n];
        }
        H[n][k + 1][k] = norm(v[k + 1][n]);
        v[k + 1][n] *= (1.0 / H[n][k + 1][k]);
        for (int i = 0; i < k; ++i) {
          Scalar giv_r = H[n][i][k];
          Scalar giv_h = H[n][i + 1][k];
          H[n][i][k] = sn[n][i] * giv_r + cs[n][i] * giv_h;
          H[n][i + 1][k] = cs[n][i] * giv_r - sn[n][i] * giv_h;
        }
        Scalar giv_r = H[n][k][k];
        Scalar giv_h = H[n][k + 1][k];
        Scalar co, si;
        double D = sqrt(giv_r * giv_r + giv_h * giv_h);
        if (giv_r > 0) {
          co = giv_h / D;
          si = giv_r / D;
          H[n][k][k] = D;
        } else {
          co = -giv_h / D;
          si = -giv_r / D;
          H[n][k][k] = -D;
        }
        H[n][k + 1][k] = 0;
        giv_r = e[n][k];
        e[n][k] = si * giv_r;
        e[n][k + 1] = co * giv_r;
        cs[n][k] = co;
        sn[n][k] = si;
      }
      for (int n = 0; n < ksize; ++n) {
        if (abs(e[n][k + 1]) < epsilon[K[n]] || iter == max_iter - 1 || k == M - 1) {
          k++; // one rhs is completed in inner loop or maximum of iterations
          for (int i = k - 1; i >= 0; --i) {
            Scalar s = e[n][i];
            for (int j = i + 1; j < k; ++j)
              s -= H[n][i][j] * y[n][j];
            y[n][i] = s / H[n][i][i];
          }
          t[K[n]] = 0;
          for (int i = 0; i < k; ++i)
            t[K[n]] += y[n][i] * v[i][n];
          k--;
          e.erase(e.begin() + n);
          y.erase(y.begin() + n);
          cs.erase(cs.begin() + n);
          sn.erase(sn.begin() + n);
          K.erase(K.begin() + n);
          H.erase(H.begin() + n);

          if (iter_min == -1) iter_min = iter + 1;
          d_min = min(d_min, abs(e[n][k + 1]));
          d_max = max(d_max, abs(e[n][k + 1]));
          for (int i = 0; i < k + 1; ++i)
            v[i].erase(n);
          --ksize;
          --n;
        }
      }
      ++iter;
      if (ksize == 0) break;
    }
    c = B * t;
    U += c;
    R -= A * c;
    dvec = R.norm();
  }
  for (int n = 0; n < size; ++n) { // set remaining rhs
    u[N[n]] = U[n];
    r[N[n]] = R[n];
  }
  if (iter == max_iter) {
    if (iter_min == -1) iter_min = iter;
    for (int n = 0; n < size; ++n) {
      d_min = min(d_min, dvec[n]);
      d_max = max(d_max, dvec[n]);
    }
  }
  d = d_max;
}