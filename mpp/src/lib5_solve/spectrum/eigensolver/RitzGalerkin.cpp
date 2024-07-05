#include "RitzGalerkin.hpp"

void RitzGalerkin::compute(
    Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    Vector &tmp) {
  mout.StartBlock("Ritz Galerkin");
  int R = u.size();

  Vector residual(u[0]); // also used as tmp vector

  SymRMatrix a(R), b(R);
  RMatrix e(R);

  int step = 0;
  int r_eps = 0;
  while (step < maxstep) {
    RitzStep(u, lambda, A, B, a, b, e, residual, tmp);
    if (r_eps == R) break;
    output(lambda, step);
    ++step;
    r_eps = 0;
    for (int r = 0; r < R; ++r) {
      double defect = ritzDefect(u[r], residual, lambda[r], A, B, tmp);
      if (defect < eps) ++r_eps;
      else u[r] -= (*solver) * residual;
    }
  }
  outputFinish(lambda, step);
  mout.EndBlock();
}
