#ifndef MLUQ_PDE_MODEL_H
#define MLUQ_PDE_MODEL_H

#include "EllipticPDESolver.hpp"
#include "InverseProblem.hpp"
#include "ObservationOperator.hpp"




class EllipticPDEModel : public EllipticPDESolver {
protected:
  std::shared_ptr<EllipticInverseProblem> problem;

  ObservationOperator observation_operator;

public:
  EllipticInverseProblem &GetProblem();

  bool run(RVector &observation);

  EllipticPDEModel(PDESolverConfig conf);

  virtual void Run(RVector &solution);
};

class TestModel : public EllipticPDEModel {
public:
  TestModel(const PDESolverConfig &conf) ;

  void Run(RVector &solution);
};

EllipticPDEModel *CreatePDEModel(PDESolverConfig conf);

std::shared_ptr<EllipticPDEModel> CreatePDEModelShared(PDESolverConfig conf);

#endif // MLUQ_PDE_MODEL_H