#ifndef BICGSTAB_HPP
#define BICGSTAB_HPP

#include "LinearSolver.hpp"

class BiCGStab : public LinearSolver {
  int restart;
public:
  explicit BiCGStab(Preconditioner *pc, const std::string &prefix = "Linear") :
      LinearSolver(pc, prefix), restart(100) {
    Config::Get(prefix + "Restart", restart);
  }

  explicit BiCGStab(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix), restart(100) {
    Config::Get(prefix + "Restart", restart);
  }

  void Restart(int &cnt, Scalar &rho, Scalar &alpha, Scalar &omega, Vector &a, Vector &b) const {
    cnt = 0;
    rho = alpha = omega = 1;
    a = 0;
    b = 0;
  }

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    Vector p(0.0, u), v(0.0, u), y(0.0, u), t(0.0, u), r0(r);
    double rho = 0.0, alpha = 0.0, omega = 0.0;
    double d = d0;
    int cnt;
    Restart(cnt, rho, alpha, omega, p, v);
    int iter = 0;
    for (; iter < max_iter; ++iter, ++cnt) {
      iteration.push_back(d);
      if (d < epsilon) break;
      if (iter % printSteps == 0) { vout(1) << iteration; }
      Scalar rho_new = r0 * r;
      if ((abs(rho) < 1e-8 * abs(rho_new) && iter) || abs(omega) < 1e-8 * abs(alpha)) {
        vout(2) << "breakdown occured (rho = " << rho << " , omega = " << omega << ")" << endl;
        Restart(cnt, rho, alpha, omega, p, v);
        continue;
      } else if (cnt >= restart) Restart(cnt, rho, alpha, omega, p, v);
      Scalar beta = (rho_new / rho) * (alpha / omega);
      p -= omega * v;
      p *= beta;
      p += r;
      y = B * p;
      v = A * y;
      Scalar h = r0 * v;
      if (abs(h) < 1e-8 * abs(rho_new)) {
        vout(2) << "breakdown occured (h=" << h << ")" << endl;
        Restart(cnt, rho, alpha, omega, p, v);
        continue;
      }
      alpha = rho_new / h;
      u += alpha * y;
      r -= alpha * v;
      d = norm(r);
      ++iter;
      if (d < epsilon) break;
      y = B * r;
      t = A * y;
      omega = (t * r) / (t * t);
      u += omega * y;
      r -= omega * t;
      rho = rho_new;
      d = norm(r);
    }
  }

  string Name() const { return "BiCGStab"; }
};

class BiCGStab2 : public LinearSolver {
public:
  explicit BiCGStab2(const std::string &prefix = "Linear") : LinearSolver(GetPC(), prefix) {}

  explicit BiCGStab2(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix) {}

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    double d = d0;
    Vector y = r;
    Scalar delta = y * r;
    Vector s = r;
    int iter = 0;
    for (; iter < max_iter; ++iter) {
      vout(2) << iter << " delta " << delta << endl;
      if (abs(delta) < epsilon) {
        y = r;
        delta = y * r;
      }
      iteration.push_back(d);
      if (d < epsilon) break;
      if (iter % printSteps == 0) { vout(1) << iteration; }
      Vector Bs = B * s;
      Vector ABs = A * Bs;
      Scalar phi = y * ABs / delta;
      vout(2) << iter << " phi " << delta << endl;
      if (abs(phi) < epsilon) break;
      Scalar omega = 1.0 / phi;
      u += omega * Bs;
      Vector w = r;
      w -= omega * ABs;
      Vector Bw = B * w;
      Vector ABw = A * Bw;
      Scalar chi = (ABw * w) / (ABw * ABw);
      vout(2) << iter << " chi " << chi << endl;
      if (abs(chi) < epsilon) break;
      r = w;
      u += chi * Bw;
      r -= chi * ABw;
      s -= chi * ABs;
      Scalar delta_new = y * r;
      Scalar psi = -omega * delta_new / (delta * chi);
      vout(2) << iter << " psi " << psi << endl;
      if (abs(psi) < epsilon) break;
      delta = delta_new;
      s *= -psi;
      s += r;
      d = norm(r);
    }
  }

  string Name() const override { return "BiCGStab2"; }
};

#endif // BICGSTAB_HPP
