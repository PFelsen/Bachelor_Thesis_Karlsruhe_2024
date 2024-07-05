#ifndef MIXEDMAIN_HPP
#define MIXEDMAIN_HPP

#include "Algebra.hpp"
#include "EllipticProblems.hpp"

class MixedEllipticAssemble;

class IProblem;

class MixedMain {
private:
  int verbose = 1;

  std::shared_ptr<IEllipticProblem> problem;

  std::shared_ptr<MixedEllipticAssemble> assemble;

  std::map<size_t, std::shared_ptr<Vector>> solutionOnLevel;

public:
  double EvaluateNormalFlux(const LevelPair& level, const Cell &c, int face) const;

  VectorField EvaluateCellFlux(const LevelPair& level, const Cell &c) const;

  void Solve();

  MixedMain(const std::string &problemName);

  void PlotVtu(const Vector &u) const;

  IProblem& GetFluxProblem();

};

#endif //MIXEDMAIN_HPP
