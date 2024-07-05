#ifndef SSOR_HPP
#define SSOR_HPP

#include "PreconditionerBase.hpp"
#include "RMatrix.hpp"

class SSOR :
    public Preconditioner { // Has to stay for historic reasons (Niethammer got famous with it)
  vector<RMatrix> D;
  int N;
  double omega;
  const Matrix *AA;
public:
  SSOR() : omega(1) { Config::Get("omega", omega); }

  SSOR(double omega) : omega(omega) {} // constructor to ensure omega=1 in the SGS case

  void Construct(const Matrix &_A) override;

  void Destruct() { D.clear(); }

  virtual ~SSOR() { Destruct(); }

  void multiply(Vector &u, const Vector &b) const override;

  string Name() const { return "SSOR"; }
};

#endif // SSOR_HPP
