#ifndef HYBRIDMAIN_HPP
#define HYBRIDMAIN_HPP

#include "Algebra.hpp"
#include "EllipticProblems.hpp"

class HybridEllipticAssemble;

class HybridMain {
private:
  int verbose = 1;

  std::shared_ptr<IEllipticProblem> problem;

  std::shared_ptr<HybridEllipticAssemble> assemble;

  std::shared_ptr<Vector> normalFlux;

  std::shared_ptr<Vector> solution;

public:
  double EvaluateNormalFlux(const Cell &c, int face) const;

  VectorField EvaluateCellFlux(const Cell &c) const;

  void Solve();

  HybridMain(const std::string &problemName);

  const Vector& NormalFlux() const { return *normalFlux; }

  const Vector& Solution() const { return *solution; }

  void PlotVtu(const Vector &u) const;

  const IProblem& GetFluxProblem() const;

};

#endif //HYBRIDMAIN_HPP
