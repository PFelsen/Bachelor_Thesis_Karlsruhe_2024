#ifndef EIGENPAIR_H
#define EIGENPAIR_H

#include "Algebra.hpp"
#include "Spectrum.hpp"

typedef Vector Eigenfct;

typedef Vectors Eigenfcts;

class Eigenpair : public VectorMatrixBase {
  Eigenfct U;
  Eigenvalue lambda;
public:
  Eigenpair(std::shared_ptr<const IDiscretization> disc, int spaceLevel = -1, int timeLevel = -1);

  Eigenpair(std::shared_ptr<const IDiscretization> disc, LevelPair levels);

  Eigenpair(const Eigenfct &);

  Eigenpair(const Eigenpair &ep);

  Eigenpair(const Eigenfct &U, const Eigenvalue &lambda);

#ifdef BUILD_IA
  void SetIADisc(std::shared_ptr<const IAIDiscretization> disc) {
    iadisc = disc;
    U.SetIADisc(disc);
  }
#endif

  Eigenpair &operator=(const Eigenpair &ep);

  Eigenpair &operator+=(const double &shift);

  Eigenpair &operator-=(const double &shift);

  const Eigenfct &getEigenfct() const { return U; }

  const Eigenvalue &getEigenvalue() const { return lambda; }

  friend inline std::ostream &operator<<(std::ostream &s, const Eigenpair &ep) {
    return s << ep.lambda;
  }

  friend inline Saver &operator<<(Saver &s, Eigenpair &ep) { return s << ep.lambda << ep.U; }

  friend inline Saver &operator<<(Saver &s, const Eigenpair &ep) { return s << ep.lambda << ep.U; }

  friend inline Loader &operator>>(Loader &l, Eigenpair &ep) { return l >> ep.lambda >> ep.U; }
};

class Eigenpairs : public VectorMatrixBase {
private:
  Eigenfcts U;
  Eigenvalues lambda;
public:
  Eigenpairs(int N, std::shared_ptr<const IDiscretization> disc, int spaceLevel = -1,
             int timeLevel = -1);

  Eigenpairs(int N, std::shared_ptr<const IDiscretization> disc, LevelPair levels);

  Eigenpairs(int N, const Eigenpairs &EP);

  Eigenpairs(const Eigenpairs &EP);

  /// selects part of EP (indices included)
  Eigenpairs(const Eigenpairs &EP, int startIdx, int endIdx = -1);

  void set(const Eigenpairs &EP);

#ifdef BUILD_IA
  void SetIADisc(std::shared_ptr<const IAIDiscretization> disc) {
    iadisc = disc;
    U.SetIADisc(disc);
  }
#endif

  Eigenpairs &operator=(const Eigenpairs &EP);

  Eigenpairs &operator+=(const double shift);

  Eigenpairs &operator-=(const double shift);

  const Eigenpair getEigenpair(int n) const { return Eigenpair(U[n], lambda[n]); }

  const Eigenfcts &getEigenfcts() const { return U; }

  const Eigenfct &LastEigenfct() const { return U.back(); }

  const Eigenvalues &getEigenvalues() const { return lambda; }

  const Eigenfct &operator()(int n) const { return U[n]; }

  const Eigenvalue &operator[](unsigned int n) const { return lambda[n]; }

  const Eigenvalue &LastEigenvalue() const { return lambda.back(); }

  int size() const { return U.size(); }

  void resize(int N);

  friend inline std::ostream &operator<<(std::ostream &s, const Eigenpairs &EP) {
    return s << EP.lambda;
  }

  friend inline Saver &operator<<(Saver &s, Eigenpairs &EP) { return s << EP.lambda << EP.U; }

  friend inline Saver &operator<<(Saver &s, const Eigenpairs &EP) { return s << EP.lambda << EP.U; }

  friend inline Loader &operator>>(Loader &l, Eigenpairs &EP) { return l >> EP.lambda >> EP.U; }

  template<bool VERIFIED>
  friend class IAEigenvalueMethod;

  template<bool VERIFIED>
  friend class IAHomotopyMethod;
};

#endif // EIGENPAIR_H
