#ifndef ITRANSPORTASSEMBLE_HPP
#define ITRANSPORTASSEMBLE_HPP

#include "TransportProblems.hpp"
#include "IDiscretization.hpp"
#include "TimeSeries.hpp"
#include "Assemble.hpp"


typedef std::pair<double, double> RatePair;

class ITransportAssemble : public ILinearTimeAssemble {
protected:
  const ITransportProblem &problem;

  int plotting = 0;

public:
  ITransportAssemble(const ITransportProblem &problem) : problem(problem) {
    Config::Get("VtuPlot", plotting);
  }

  const ITransportProblem &GetProblem() const { return problem; }

  virtual const IDiscretization &GetDisc() const = 0;

  virtual std::shared_ptr<const IDiscretization> GetSharedDisc() const = 0;

  void SetTimeSeries(const Vector &u) {
    ResetTime(problem.GetStartTime(), problem.GetEndTime(),
              problem.GetStepSize(u.GetMesh().Level().space));
  }

  virtual double Energy(const Vector &u) const = 0;

  virtual double Mass(const Vector &u) const = 0;

  virtual double Error(const Vector &u) const = 0;

  virtual RatePair InflowOutflow() const = 0;
  
  virtual RatePair InflowOutflow(const Vector &u) const = 0;

  virtual double MaxFlux(const Vector &u) const = 0;

  virtual double Flow(const Vector &u) const = 0;
};

struct PDESolverConfig;

ITransportAssemble *
CreateTransportAssemble(const ITransportProblem &problem, const PDESolverConfig &conf);

std::unique_ptr<ITransportAssemble>
CreateTransportAssembleUnique(const ITransportProblem &problem, const PDESolverConfig &conf);

std::shared_ptr<ITransportAssemble>
CreateTransportAssembleShared(const ITransportProblem &problem, const PDESolverConfig &conf);

#endif //ITRANSPORTASSEMBLE_HPP
