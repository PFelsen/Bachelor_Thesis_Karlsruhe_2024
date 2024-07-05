

#ifndef CG_HPP
#define CG_HPP

#include <complex>
#include "LinearSolver.hpp"

class CG : public LinearSolver {
public:
  CG(Preconditioner *pc, const std::string &prefix = "Linear") : LinearSolver(pc, prefix) {}

  CG(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix) {}

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    double d = d0;
    Vector c(u);
    Vector p(u);
    Vector t(u);
    p = 0;
    int iter = 0;
    Scalar rho_0 = 1;
    for (; iter < max_iter; ++iter) {
      if (d < epsilon) break;
      vout(1) << iteration;
      vout(10) << "r " << r << endl;
      c = B * r;
      vout(10) << "c " << c << endl;
      Scalar rho_1 = r * c;
      vout(2) << "rho " << r * c << " " << c * r << endl;
      p *= (rho_1 / rho_0);
      p += c;
      rho_0 = rho_1;
      t = A * p;
      Scalar alpha = rho_0 / (p * t);
      vout(3) << "alpha " << p * t << " " << t * p << endl;
      u += alpha * p;
      vout(10) << "u " << u << endl;
      r -= alpha * t;
      d = r.norm();
      iteration.push_back(d);
    }
  }

  string Name() const override { return "PCG"; }
};

class CGNE : public LinearSolver {
public:
  CGNE(Preconditioner *pc, const std::string &prefix = "Linear") : LinearSolver(pc, prefix) {}

  CGNE(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix) {}

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    Vector Atr(r);
    Atr = r * A;
    double d = Atr.norm();
    double d_0 = d;
    double eps = defaultEpsilon + defaultReduction * d;
    Vector Bc(u);
    Vector c(u);
    Vector p(u);
    Vector Ap(u);
    Vector t(u);
    p = 0;
    Scalar rho_0 = 1;
    int iter = 0;
    for (; iter < max_iter; ++iter) {
      if (d < eps) break;
      vout(1) << iteration;
      vout(10) << "r " << Atr << endl;
      Bc = B * Atr;
      Bc.Collect();
      c = Bc * B;
      vout(10) << "c " << c << endl;
      Scalar rho_1 = Atr * c;
      vout(2) << "rho " << Atr * c << " " << c * Atr << endl;
      p *= (rho_1 / rho_0);
      p += c;
      rho_0 = rho_1;
      Ap = A * p;
      Ap.Accumulate();
      t = Ap * A;
      Scalar alpha = rho_0 / (p * t);
      vout(3) << "alpha " << p * t << " " << t * p << endl;
      u += alpha * p;
      vout(10) << "u " << u << endl;
      Atr -= alpha * t;
      d = Atr.norm();
      iteration.push_back(d);
    }
  }

  string Name() const { return "PCGNE"; }
};

// TODO: Replace ScalarMatrix and use Spectrum for SpectralBounds!
#include "MINRES.hpp"


extern "C" void zhegv_(int *ITYPE, char *JOBZ, char *UPLO, int *N, void *A, int *LDA, void *B,
                       int *LDB, void *W, void *WORK, int *LWORK, void *RWORK, int *INFO);

std::pair<double, double> SpectralBounds(int n, const RMatrix &a) {
  int verbose = -1;
  Config::Get("SpectralBoundsVerbose", verbose, true);
  if (n == 1) {
    if (a[0][0] == 0.0) return std::pair<double, double>(1, 1);
    return std::pair<double, double>(a[0][0], a[0][0]);
  }
  char jobz = 'N';
  char uplo = 'U';
  int itype = 1;
  int lwork = 2 * n;
  int info = -1;
  double *w = new double[n];
  std::complex<double> *work = new std::complex<double>[2 * n];
  double *rwork = new double[3 * n];
  std::complex<double> *A = new std::complex<double>[n * n];
  std::complex<double> *B = new std::complex<double>[n * n];
  vout(10) << "A: " << DOUT(n) << endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i * n + j] = a[i][j];
      vout(10) << a[i][j] << ", ";
      B[i * n + j] = (i == j);
    }
    vout(10) << endl;
  }
  zhegv_(&itype, &jobz, &uplo, &n, A, &n, B, &n, w, work, &lwork, rwork, &info);

  vout(2) << "SpectralBounds INFO: " << info << endl;

  double w_min = w[0];
  double w_max = w[n - 1];
  delete[] w;
  delete[] work;
  delete[] rwork;
  delete[] A;
  delete[] B;

  if (info != 0) { Warning("Error in Lapack-Routine to calc Eigenvalues.") return {1, -1}; }
  return std::pair<double, double>(w_min, w_max);
}

class CGX : public LinearSolver {
public:
  CGX(Preconditioner *pc, const std::string &prefix = "Linear") : LinearSolver(pc, prefix) {}

  CGX(std::unique_ptr<Preconditioner> pc, const std::string &prefix = "Linear") :
      LinearSolver(std::move(pc), prefix) {}

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const override {
    double d = r.norm();
    double eps = defaultEpsilon + defaultReduction * d;
    Vector c(0.0, u);
    Vector p(0.0, u);
    Vector t(0.0, u);
    Scalar rho_0 = 1;
    RMatrix H(0.0, max_iter + 1);
    Vectors U(max_iter + 1, u);
    Vectors AU(max_iter + 1, r);
    U[0] = B * r;
    AU[0] = A * U[0];
    Scalar h_10 = sqrt(U[0] * AU[0]);
    U[0] *= 1.0 / h_10;
    AU[0] *= 1.0 / h_10;
    std::pair<double, double> Kappa;
    int iter = 0;
    for (; iter < max_iter; ++iter) {
      Vector BAu = B * AU[iter];
      U[iter + 1] = BAu;
      for (int i = 0; i <= iter; ++i) {
        H[i][iter] = BAu * AU[i];
        U[iter + 1] -= H[i][iter] * U[i];
      }
      AU[iter + 1] = A * U[iter + 1];
      H[iter + 1][iter] = sqrt(AU[iter + 1] * U[iter + 1]);
      U[iter + 1] *= 1.0 / H[iter + 1][iter];
      AU[iter + 1] *= 1.0 / H[iter + 1][iter];
      Kappa = SpectralBounds(iter + 1, H);
      if (d < eps) break;
      vout(1) << "d(" << iter << ")= " << d << " kappa=" << Kappa.second / Kappa.first << " spec=["
              << Kappa.first << "," << Kappa.second << "]" << endl;
      vout(10) << "r " << r << endl;
      c = B * r;
      vout(10) << "c " << c << endl;
      Scalar rho_1 = r * c;
      vout(2) << "rho " << r * c << " " << c * r << endl;
      p *= (rho_1 / rho_0);
      p += c;
      rho_0 = rho_1;
      t = A * p;
      Scalar alpha = rho_0 / (p * t);
      vout(3) << "alpha " << p * t << " " << t * p << endl;
      u += alpha * p;
      vout(10) << "u " << u << endl;
      r -= alpha * t;
      d = r.norm();
      iteration.push_back(d);
    }
    double s = 0;
    for (int i = 1; i < iter; ++i) {
      for (int j = i - 1; j < i; ++j) {
        s += abs(H[i][j] - H[j][i]);
      }
    }
    vout(1) << "kappa = " << Kappa.second / Kappa.first << " spec=[" << Kappa.first << ","
            << Kappa.second << "]"
            << " |skew H| = " << s << endl;
    vout(2) << " H = " << endl << H << endl;
  }

  string Name() const { return "PCGX"; }
};

#endif // CG_HPP