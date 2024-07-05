//
// Created by lstengel on 21.07.22.
//

#ifndef TUTORIAL_VECTORVALUEDMAIN_HPP
#define TUTORIAL_VECTORVALUEDMAIN_HPP

#include <utility>

#include "IVectorValuedAssemble.hpp"
#include "EllipticProblems.hpp"
#include "Newton.hpp"
#include "PDESolver.hpp"

class VectorValuedMain : public PDESolver<IVectorValuedProblem> {
private:
  std::shared_ptr<Newton> newton; // Todo replace with linear solver

  std::shared_ptr<IVectorValuedAssemble> assemble;
protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<IVectorValuedProblem> problem) override;

public:
  explicit VectorValuedMain(const PDESolverConfig &conf) : GenericPDESolver(conf),
  newton(CreateNewton(conf.linearSolver, conf.preconditioner))  {

  }

  void PrintValues(const Solution &solution) override;

  double EvaluateQuantity(const Vector &u, const string &quantity) const;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};


#endif //TUTORIAL_VECTORVALUEDMAIN_HPP
