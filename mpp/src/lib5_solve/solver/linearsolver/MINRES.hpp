#ifndef MINRES_HPP
#define MINRES_HPP

#include "LinearSolver.hpp"

class MINRES : public LinearSolver {
  int M;
public:
  MINRES(Preconditioner *pc, const std::string &prefix = "Linear") : LinearSolver(pc, prefix) {
    M = 5;
  }

  MINRES(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix) {
    M = 5;
  }

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    double d = d0;
    RMatrix H(M + 1, M);
    RVector e(M);
    Vectors v(M + 1, u);
    Vector c(r), t(r), y(r);
    int iter = 0;
    for (; iter < max_iter;) {
      if (d < epsilon) break;
      vout(1) << iteration;
      v[0] = (1 / d) * r;
      int k = 0;
      for (; k < M; ++k) {
        c = B * v[k];
        t = A * c;
        v[k + 1] = t;
        for (int i = 0; i <= k; ++i) {
          H[i][k] = t * v[i];
          v[k + 1] -= H[i][k] * v[i];
        }
        H[k + 1][k] = norm(v[k + 1]);
        vout(3) << "GMRES: " << k << " |v_k+1| " << H[k + 1][k] << endl;
        v[k + 1] *= (1.0 / H[k + 1][k]);
      }
      for (int i = k - 1; i >= 0; --i) {
        Scalar s = e[i];
        for (int j = i + 1; j < k; ++j)
          s -= H[i][j] * y[j];
        y[i] = s / H[i][i];
      }
      t = 0;
      for (int i = 0; i < k; ++i)
        t += y[i] * v[i];
      c = B * t;
      u += c;
      r -= A * c;
      d = norm(r);
      iteration.push_back(d);
    }
  }

  string Name() const { return "MINRES"; }
};

#endif // MINRES_HPP
