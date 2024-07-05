#include "Eigenpair.hpp"

Eigenpair::Eigenpair(std::shared_ptr<const IDiscretization> disc, int spaceLevel, int timeLevel) :
    Eigenpair(disc, {spaceLevel, timeLevel}) {}

Eigenpair::Eigenpair(std::shared_ptr<const IDiscretization> disc, LevelPair levels) :
    VectorMatrixBase(disc, levels), U(disc, levels), lambda(0) {}

Eigenpair::Eigenpair(const Eigenfct &U) : VectorMatrixBase(U), U(U), lambda(0) {}

Eigenpair::Eigenpair(const Eigenpair &ep) : VectorMatrixBase(ep), U(ep.U), lambda(ep.lambda) {}

Eigenpair::Eigenpair(const Eigenfct &U, const Eigenvalue &lambda) :
    VectorMatrixBase(U), U(U), lambda(lambda) {}

Eigenpair &Eigenpair::operator=(const Eigenpair &ep) {
  U = ep.U;
  lambda = ep.lambda;
#ifdef BUILD_IA
  iadisc = ep.iadisc;
#endif
  return *this;
}

Eigenpair &Eigenpair::operator+=(const double &shift) {
  lambda += shift;
  return *this;
}

Eigenpair &Eigenpair::operator-=(const double &shift) {
  lambda -= shift;
  return *this;
}

Eigenpairs::Eigenpairs(int N, std::shared_ptr<const IDiscretization> disc, int spaceLevel,
                       int timeLevel) : Eigenpairs(N, disc, {spaceLevel, timeLevel}) {}

Eigenpairs::Eigenpairs(int N, std::shared_ptr<const IDiscretization> disc, LevelPair levels) :
    VectorMatrixBase(disc, levels), U(N, disc, levels), lambda(N) {}

Eigenpairs::Eigenpairs(int N, const Eigenpairs &EP) : VectorMatrixBase(EP), U(N, EP.U), lambda(N) {}

Eigenpairs::Eigenpairs(const Eigenpairs &EP) : VectorMatrixBase(EP), U(EP.U), lambda(EP.lambda) {}

Eigenpairs::Eigenpairs(const Eigenpairs &EP, int startIdx, int endIdx) :
    VectorMatrixBase(EP), U(EP.U, startIdx, endIdx),
    lambda(EP.lambda, startIdx, startIdx + U.size()) {}

void Eigenpairs::set(const Eigenpairs &EP) {
  for (int i = 0; i < min(size(), EP.size()); ++i) {
    U[i] = EP.U[i];
    lambda[i] = EP.lambda[i];
  }
  for (int i = min(size(), EP.size()); i < size(); ++i) {
    U[i] = 0;
    lambda[i] = 0;
  }
#ifdef BUILD_IA
  iadisc = EP.iadisc;
#endif
}

Eigenpairs &Eigenpairs::operator=(const Eigenpairs &EP) {
  U.resize(EP.size());
  lambda.resize(EP.size());
  U = EP.U;
  lambda = EP.lambda;
#ifdef BUILD_IA
  iadisc = EP.iadisc;
#endif
  return *this;
}

Eigenpairs &Eigenpairs::operator+=(const double shift) {
  lambda += shift;
  return *this;
}

Eigenpairs &Eigenpairs::operator-=(const double shift) {
  lambda -= shift;
  return *this;
}

void Eigenpairs::resize(int N) {
  if (N == U.size()) return;
  U.resizeData(N);
  lambda.resize(N);
}
