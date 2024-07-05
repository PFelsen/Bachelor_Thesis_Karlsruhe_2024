#ifndef MLUQ_STOCHASTIC_INVERSE_PROBLEM_H
#define MLUQ_STOCHASTIC_INVERSE_PROBLEM_H

#include "IEllipticAssemble.hpp"
#include "IProblem.hpp"
#include <tuple>
#include "KLExpansion.hpp"
#include "ProposalGenerators.hpp"


class ProposalInput {
public:
  ProposalInput(ProposalConfig conf) {}

  virtual KLExpansion &GetRandomField() = 0;
};

class ProposalInput0D : public ProposalInput {
public:
  std::unique_ptr<KLExpansion> random_field_generator;

  ProposalInput0D(ProposalConfig conf) : ProposalInput(conf) {
    random_field_generator = CreateRandomFieldUnique(conf.field_name);
  }

  double EvaluateSample(const Point &x) const { return random_field_generator->EvalSample(x); }

  KLExpansion &GetRandomField() override { return *random_field_generator; }
};

class ProposalInput2D : public ProposalInput {
private:
  std::unique_ptr<KLExpansion> random_field_generator;
public:
  ProposalInput2D(ProposalConfig conf) : ProposalInput(conf) {
    random_field_generator = CreateRandomFieldUnique(conf.field_name);
  }

  double EvaluateSample(const Point &x) const { return random_field_generator->EvalSample(x); }

  KLExpansion &GetRandomField() override { return *random_field_generator; }
};

class EllipticInverseProblem : public IEllipticProblem {
public:
  EllipticInverseProblem() : IProblem("Interval") {}

  virtual ProposalInput &GetInput() = 0;

  Scalar Load(const Point &x, const cell &c) const override { return 0.0; }

  Scalar Solution(const Point &x) const override { return Permeability(x)[0][0]; }

  VectorField Flux(const Point &x) const override { return {-1.0, 0.0}; }

  Tensor Permeability(const Point &x) const override { return One; }

  std::string Name() const override { return "Stupid virtual workaround (temporary)"; }
};

class TestWithLaplace1D : public EllipticInverseProblem {
public:
  ProposalInput0D input;

  ProposalInput &GetInput() override { return input; }

  std::string Name() const override { return "Test Problem (Laplace 1D)"; }

  TestWithLaplace1D(ProposalConfig conf) : input(conf), IProblem("Interval") {}

  bool HasExactSolution() const override { return true; }

  Scalar Load(const Point &x, const cell &c) const override { return 0.0; }

  Scalar Solution(const Point &x) const override { return Permeability(x)[0][0]; }

  VectorField Flux(const Point &x) const override { return {-1.0, 0.0}; }

  Tensor Permeability(const Point &x) const override { return input.EvaluateSample(x) * One; }
};

class TestWithLaplace2D : public EllipticInverseProblem {
public:
  std::string Name() const override { return "Test Problem (Laplace 2D)"; }

  ProposalInput0D input;

  ProposalInput &GetInput() override { return input; }

  TestWithLaplace2D(ProposalConfig conf) : input(conf), IProblem("UnitSquare") {}

  bool HasExactSolution() const override { return true; }

  Scalar Load(const Point &x, const cell &c) const override { return 0.0; }

  Scalar Solution(const Point &x) const override { return -x[1]; }

  VectorField Flux(const Point &x) const override { return {0.0, -1.0}; }

  Tensor Permeability(const Point &x) const override { return input.EvaluateSample(x) * One; }
};

class RainLaplace2D : public EllipticInverseProblem {
public:
  std::string Name() const override { return "Rain subsurface diffusion problem (Laplace 2D)"; }

  ProposalInput2D input;

  ProposalInput &GetInput() override { return input; }

  RainLaplace2D(ProposalConfig conf) : input(conf), IProblem("UnitSquare") {}

  bool HasExactSolution() const override { return true; }

  Scalar Load(const Point &x, const cell &c) const override { return 0.0; }

  Scalar Solution(const Point &x) const override { return -x[1]; }

  VectorField Flux(const Point &x) const override { return {0.0, -1.0}; }

  Tensor Permeability(const Point &x) const override { return input.EvaluateSample(x) * One; }
};

EllipticInverseProblem *CreateEllipticInverse(std::string problem, ProposalConfig conf);

std::shared_ptr<EllipticInverseProblem> CreateEllipticInverseShared(std::string problem,
                                                                    ProposalConfig conf);

#endif // MLUQ_STOCHASTIC_INVERSE_PROBLEM_H