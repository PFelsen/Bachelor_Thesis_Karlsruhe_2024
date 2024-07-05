#ifndef IACOUSTICASSEMBLE_HPP
#define IACOUSTICASSEMBLE_HPP

#include "AcousticProblems.hpp"
#include "Assemble.hpp"

struct AcousticErrors {
  double l1 = 0.0;

  double l2 = 0.0;

  double inf = 0.0;
};

struct AcousticNorms {
  double l1 = 0.0;

  double l2 = 0.0;

  double energy = 0.0;
};

typedef std::pair<double, double> MinMaxP;

class IAcousticAssemble : public ILinearTimeAssemble {
protected:
  const AcousticProblem &problem;

  int plotting = 0;

public:
  IAcousticAssemble(const AcousticProblem &problem) : problem(problem) {
    Config::Get("VtuPlot", plotting);
  }

  void SetTimeSeries(const Vector &u) {
    ResetTime(problem.GetStartTime(), problem.GetEndTime(),
              problem.GetStepSize(u.GetMesh().Level().space));
  }

  const AcousticProblem &GetProblem() const { return problem; }

  virtual const IDiscretization &GetDisc() const = 0;

  virtual std::shared_ptr<const IDiscretization> GetSharedDisc() const = 0;

  virtual const char *Name() const = 0;

  virtual double Energy(const Vector &u) const = 0;

  virtual AcousticNorms ComputeNorms(const Vector &u) const = 0;

  virtual AcousticErrors ComputeErrors(const Vector &u) const = 0;

  virtual MinMaxP MinMaxPressure(const Vector &u) const = 0;
};

struct PDESolverConfig;

IAcousticAssemble *
CreateAcousticAssemble(const AcousticProblem &problem, const PDESolverConfig &conf);

std::unique_ptr<IAcousticAssemble>
CreateAcousticAssembleUnique(const AcousticProblem &problem, const PDESolverConfig &conf);

std::shared_ptr<IAcousticAssemble>
CreateAcousticAssembleShared(const AcousticProblem &problem, const PDESolverConfig &conf);

#endif //IACOUSTICASSEMBLE_HPP
