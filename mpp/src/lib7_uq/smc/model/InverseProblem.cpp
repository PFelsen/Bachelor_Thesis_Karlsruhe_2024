#include "InverseProblem.hpp"


EllipticInverseProblem *CreateEllipticInverse(std::string problem, ProposalConfig conf) {
  if (problem == "Test2D") { return new TestWithLaplace2D(conf); }
  if (problem == "Test1D") { return new TestWithLaplace1D(conf); }
  if (problem == "RainLaplace2D") { return new RainLaplace2D(conf); }
  throw std::invalid_argument("Inverse problem type " + problem + " cannot be found.");
}

std::shared_ptr<EllipticInverseProblem> CreateEllipticInverseShared(std::string problem,
                                                                    ProposalConfig conf) {
  return std::shared_ptr<EllipticInverseProblem>(CreateEllipticInverse(problem, conf));
}