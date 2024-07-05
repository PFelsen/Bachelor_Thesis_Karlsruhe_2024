#ifndef REACTIONPDESOLVER_HPP
#define REACTIONPDESOLVER_HPP

#include "ReactionProblems.hpp"
#include "DGReactionAssemble.hpp"
#include "PGReactionAssemble.hpp"

#include "PDESolver.hpp"
#include "Newton.hpp"

class ReactionPDESolver : public PDESolver<IReactionProblem> {
private:
  std::shared_ptr<IReactionAssemble> assemble;
protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<IReactionProblem> problem) override;

public:
  explicit ReactionPDESolver(const PDESolverConfig &conf) : GenericPDESolver(conf) {
  }

  std::string Name() const override { return "ParabolicPDESolver"; }

  std::map<std::string, double> ComputeValues(const Vector &u) const override;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};

#endif //REACTIONPDESOLVER_HPP
