#ifndef TUTORIAL_FGMRES_H
#define TUTORIAL_FGMRES_H

#include "LinearSolver.hpp"

class FGMRES : public LinearSolver {
  int M = 100;
public:
  FGMRES(Preconditioner *pc, const std::string &prefix = "Linear") : LinearSolver(pc, prefix) {
    Config::Get("restart", M);
  }

  FGMRES(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix) {
    Config::Get("restart", M);
  }

  string Name() const { return "FGMRES"; }

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    double d = d0;
    Scalar H[M + 1][M], e[M + 1], y[M], cs[M], sn[M];
    LazyVectors v(M + 1, u);
    LazyVectors c(M + 1, u);
    Vector t(r);
    int iter = 0;
    for (; iter < max_iter;) {
      if ((iter >= min_iter) && (d < epsilon)) break;
      vout(1) << "d(" << iter << ")= " << d << "\n";
      v[0] = (1 / d) * r;
      e[0] = d;
      for (int i = 1; i < M; ++i)
        e[i] = 0;
      int k = 0;
      for (; (k < M) && (iter < max_iter); ++k) {
        c[k] = B * v[k];
        t = A * c[k];
        vout(2) << k << " c = B v " << norm(c[k]) << " Ac " << norm(t) //  << " with " << B.Name()
                << endl;
        v[k + 1] = t;
        for (int i = 0; i <= k; ++i) {
          H[i][k] = t * v[i];
          v[k + 1] -= H[i][k] * v[i];
        }
        H[k + 1][k] = norm(v[k + 1]);
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
        d = sqrt(giv_r * giv_r + giv_h * giv_h);
        vout(3) << k << " d " << d << "\n";
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

        {
          bool isAnyActive = false;
          for (auto &cb : callbacks) {
            isAnyActive |= cb->isActive(k, iter, abs(e[k + 1]));
          };
          if (isAnyActive) {
            for (int i = k - 1; i >= 0; --i) {
              Scalar s = e[i];
              for (int j = i + 1; j < k; ++j)
                s -= H[i][j] * y[j];
              y[i] = s / H[i][i];
            }
            t = 0;
            for (int i = 0; i < k; ++i) {
              t += y[i] * c[i];
            }
            Vector tt(t);
            tt += u;
            bool isAnyFinished = false;
            for (auto &cb : callbacks) {
              if (cb->isActive(k, iter, abs(e[k + 1]))) {
                double error = cb->checkCurrentSolution((const Matrix &)A, tt);
                if (error != -1) { vout(1) << "E(" << iter << ")= " << error << "\n"; }
              }
              isAnyFinished |= cb->isFinished(abs(e[k + 1]));
            }
            if (isAnyFinished) {
              u += t;
              r -= A * t;
              d = norm(r);
              return;
            }
          }
        }


        if (iter % printSteps == 0) { vout(1) << iteration; }
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
      t = 0;
      for (int i = 0; i < k; ++i)
        t += y[i] * c[i];
      u += t;
      r -= A * t;
      d = norm(r);
    }
  }
};

#endif // TUTORIAL_FGMRES_H
