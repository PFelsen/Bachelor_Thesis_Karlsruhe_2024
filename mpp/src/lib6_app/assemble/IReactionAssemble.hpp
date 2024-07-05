#ifndef IREACTIONASSEMBLE_HPP
#define IREACTIONASSEMBLE_HPP

#include "ReactionProblems.hpp"
#include "Assemble.hpp"
#include <iomanip>


typedef std::pair<double, double> RatePair;

class IReactionAssemble : public INonLinearTimeAssemble {
protected:
  const IReactionProblem &problem;

  int plotting = 0;

public:
  explicit IReactionAssemble(const IReactionProblem &problem) :
      problem(problem) {
    Config::Get("VtuPlot", plotting);
  }

  const IReactionProblem &GetProblem() const { return problem; }

  virtual const IDiscretization &GetDisc() const = 0;

  virtual std::shared_ptr<const IDiscretization> GetSharedDisc() const = 0;

  void SetTimeSeries(const Vector &u) {
    ResetTime(problem.GetStartTime(), problem.GetEndTime(),
              problem.GetStepSize(u.GetMesh().MaxMeshWidth()));
  }

  double Energy(const Vector &u) const override { return 0.0; };

  virtual double Mass(const Vector &u) const = 0;

  virtual RatePair InflowOutflow(const Vector &u) const = 0;
};

struct PDESolverConfig;

IReactionAssemble *
CreateReactionAssemble(const IReactionProblem &problem, const PDESolverConfig &conf);

std::unique_ptr<IReactionAssemble>
CreateReactionAssembleUnique(const IReactionProblem &problem, const PDESolverConfig &conf);

std::shared_ptr<IReactionAssemble>
CreateReactionAssembleShared(const IReactionProblem &problem, const PDESolverConfig &conf);

#endif //IREACTIONASSEMBLE_HPP
