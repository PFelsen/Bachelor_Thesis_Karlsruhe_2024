#ifndef GMRES2_H
#define GMRES2_H

#include "LinearSolver.hpp"

/*
 * Warning this class has no test coverage and
 * no documentation about it's functionality.
 *
 */
// Left preconditioned Version of GMRES. TODO: Rename GMRES methods
class GMRES2 : public LinearSolver {
private:
  int M = 300;
public:
  GMRES2() : LinearSolver(GetPC()) { Config::Get("restart", M); }

  void Solve(Vector &u, const Operator &A, const Operator &B, Vector &r) {
    Scalar H[M + 1][M], e[M + 1], y[M], cs[M], sn[M];
    Vectors v(M + 1, u);
    Vector c(r), t(r);
    double d_0 = norm(r);
    double d = d_0;
    double epsilon = max(config_eps, config_reduction * d);
    int iter = 0;
    for (; iter < max_iter;) {
      if ((iter >= min_iter) && (d < epsilon)) break;
      t = B * r;
      v[0] = t;
      d = sqrt(ConsistentScalarProduct(t, t));
      vout(1) << "d(" << iter << ")= " << d << endl;
      if (iter == 0) { epsilon = max(config_eps, config_reduction * d); }
      v[0] *= (1 / d);
      e[0] = d;
      for (int i = 1; i < M; ++i)
        e[i] = 0;
      int k = 0;
      for (; (k < M) && (iter < max_iter); ++k) {
        c = A * v[k];
        t = B * c;
        vout(2) << k << " c " << norm(c) << " Ac " << norm(t) << endl;
        for (int i = 0; i <= k; ++i) {
          H[i][k] = ConsistentScalarProduct(t, v[i]);
          v[k + 1] -= H[i][k] * v[i];
        }
        H[k + 1][k] = sqrt(ConsistentScalarProduct(v[k + 1], v[k + 1]));
        vout(3) << k << " |v_k+1| " << H[k + 1][k] << "\n";
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
        d = sqrt(std::real(std::conj(giv_r) * giv_r + giv_h * giv_h));
        vout(3) << k << " d " << d << "\n";
        if (std::real(giv_r) > 0) {
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
        if (iter % printSteps == 0) { vout(1) << "d(" << iter << ")= " << abs(e[k + 1]) << endl; }
        if ((iter >= min_iter) && (abs(e[k + 1]) < epsilon)) {
          k++;
          break;
        }
      }
      for (int i = k - 1; i >= 0; --i) {
        Scalar s = e[i];
        for (int j = i + 1; j < k; ++j)
          s -= H[i][j] * y[j];
        y[i] = s / H[i][i];
      }
      c = 0;
      for (int i = 0; i < k; ++i)
        c += y[i] * v[i];
      u += c;
      r -= A * c;
      d = norm(r);
    }
  }

  void Solve(Vectors &u, const Operator &A, const Operator &B, Vectors &r) {}

  string Name() const { return "GMRES2"; }
};

#endif // GMRES2_H