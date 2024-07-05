#ifndef GMRES3_H
#define GMRES3_H

#include "LinearSolver.hpp"

/*
 * Warning this class has no test coverage and
 * no documentation about it's functionality.
 *
 */
// Right preconditioned Version of with Weightmatrix, so the norm is in the W norm and
// orthogonalization with regards to W GMRES. TODO: Rename GMRES methods, unify and factorize GMRES
// => every GMRES has its own Givens rotation

class GMRES3 : public LinearSolver {
private:
  int M = 300;
public:
  GMRES3() : LinearSolver(GetPC()) { Config::Get("restart", M); }

  void Solve(Vector &u, const Operator &A, const Operator &B, Vector &r) {
    THROW("not implemented")
  }

  void Solve(Vector &u, const Operator &A, const Operator &B, const Operator &W, Vector &r) {
    Scalar H[M + 1][M], e[M + 1], y[M], cs[M], sn[M];
    LazyVectors v(M + 1, u);
    Vector c(r), t(r);
    double d_0 = norm(r);
    double d = d_0;
    double epsilon = max(config_eps, config_reduction * d);
    int iter = 0;
    for (; iter < max_iter;) {
      t = B * r;
      v[0] = t;
      d = sqrt(t * (W * t));
      if (iter == 0) { epsilon = max(config_eps, config_reduction * d); }
      vout(1) << "d(" << iter << ")= " << d << " euk " << norm(r) << endl;
      if ((iter >= min_iter) && (d < epsilon)) break;

      v[0] *= (1 / d);
      e[0] = d;
      for (int i = 1; i < M; ++i)
        e[i] = 0;
      int k = 0;
      for (; (k < M) && (iter < max_iter); ++k) {
        c = A * v[k];
        t = B * c;
        vout(2) << k << " c " << norm(c) << " Ac " << norm(t) << "\n";
        v[k + 1] = t;
        for (int i = 0; i <= k; ++i) {
          H[i][k] = (W * t) * v[i];
          v[k + 1] -= H[i][k] * v[i];
        }
        t = v[k + 1];
        H[k + 1][k] = sqrt(t * W * t);
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
    if (d < epsilon) {
      Vector Br{B * r};
      double d2 = sqrt(ConsistentScalarProduct(B * r, W * Br));
      vout(0) << "d(" << iter << ")= " << d2 << " euk " << d << " rate " << rate(iter, d_0, d)
              << endl;
      mout.EndBlock(verbose, 0);
    } else vout(-1) << "d(" << iter << ")= " << d << " rate " << rate(iter, d_0, d) << endl;
    ;
  }

  void Solve(Vectors &u, const Operator &A, const Operator &B, Vectors &r){THROW("not implemented")}

  string Name() const {
    return "GMRES3";
  }
};

#endif // TUTORIAL_GMRES3_H
