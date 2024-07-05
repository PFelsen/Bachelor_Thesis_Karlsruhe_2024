#ifndef _ASSEMBLE_H_
#define _ASSEMBLE_H_

#include "Algebra.hpp"
#include "Parallel.hpp"
#include "Plotting.hpp"
#include "TimeSeries.hpp"
#include "ctools.hpp"

class IAssemble {
protected:
  int verbose = 1;
public:
  IAssemble() { Config::Get("AssembleVerbose", verbose); }

  virtual ~IAssemble() = default;

  virtual const char *Name() const = 0;

  virtual void Initialize(Vector &u) const = 0;

  virtual double Energy(const Vector &u) const;
  ;

  virtual void Energy(const cell &c, const Vector &u,
                      double &energy) const {THROW("Virtual function not implemented")};

  virtual double Residual(const Vector &u, Vector &defect) const;

  virtual void Residual(const cell &c, const Vector &u,
                        Vector &defect) const {THROW("Virtual function not implemented")};

  virtual void Jacobi(const Vector &u, Matrix &jacobi) const;

  virtual void Jacobi(const cell &c, const Vector &u,
                      Matrix &jacobi) const {THROW("Virtual function not implemented")};

  virtual void PrintInfo() const {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 0));
  }

  virtual void PrintInfo(const Vector &u) const {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 0),
                   PrintInfoEntry("Discretization", u.GetDisc().DiscName(), 1),
                   PrintInfoEntry("Level", u.GetMesh().Level().str(), 1));
  }
};

class ILinearAssemble {
protected:
  int verbose = 1;
public:
  ILinearAssemble() { Config::Get("AssembleVerbose", verbose); }

  virtual ~ILinearAssemble() = default;

  virtual const char *Name() const = 0;

  virtual void Initialize(Vector &u) const {
    THROW("Initialize(Vector &u) not implemented in ILinearAssemble")
  }

  virtual void Initialize(Vectors &u) const {
    THROW("Initialize(Vectors &u) not implemented in ILinearAssemble")
  }

  virtual void AssembleSystem(Matrix &systemMatrix, Vector &rhs) const;

  virtual void AssembleSystem(Matrix &systemMatrix, Vectors &rhs) const;

  virtual void AssembleSystem(const cell &c, Matrix &systemMatrix, Vector &rhs) const {
    THROW("AssembleSystem(const cell &c, Matrix &systemMatrix, Vector &rhs) "
          "not implemented in ILinearAssemble")
  }

  virtual void AssembleSystem(const cell &c, Matrix &systemMatrix, Vectors &rhs) const {
    THROW("AssembleSystem(const cell &c, Matrix &systemMatrix, Vectors &rhs) "
          "not implemented in ILinearAssemble")
  }

  virtual void RHS(double t, Vector &rhs) const {
      THROW("RHS(double t, Vector &rhs) not implemented in ILinearAssemble")};

  virtual void MassMatrix(Matrix &massMatrix) const {
    THROW("MassMatrix(Matrix &massMatrix) not implemented in ILinearAssemble")
  }

  virtual void SystemMatrix(Matrix &systemMatrix) const {
    THROW("SystemMatrix(Matrix &systemMatrix) not implemented in ILinearAssemble")
  }

  virtual void PrintInfo() const {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 1));
  }

  virtual void PrintInfo(const Vector &u) const {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Discretization", u.GetDisc().DiscName(), 2),
                   PrintInfoEntry("Level", u.GetMesh().Level().str(), 2));
  }

  virtual void PrintInfo(const Vectors &u) const {
    if (u.size() < 1) {
      PrintInfo();
    } else {
      PrintInfo(u[0]);
    }
  }
};

/*
 * Todos:
 *  * Unify TimeAssemble classes by creating another Interface dealing with the
 * TimeSeries
 *  * E.g. Derive ILinearTimeAssemble from ILinearAssemble and ITimeAssemble
 *  * Rename IAssemble to INonLinearAssemble
 *  * Each Assemble Interface has to be dedicated to a specific class of solvers
 *  * Each Solver should provide a Method accepting these Interfaces
 */

class TimeInterface {
  TimeSeries timeSeries{};
protected:
  TimeInterface(TimeSeries *tSeries) : timeSeries(*tSeries) {
    delete tSeries;
    tSeries = &timeSeries;
  }
public:
  TimeInterface() {}

  explicit TimeInterface(double startTime, double endTime, double dt) :
      timeSeries(startTime, endTime, dt) {}

  explicit TimeInterface(double startTime, double endTime, int steps) :
      timeSeries(startTime, endTime, steps) {}

  void ResetTime(double startTime, double endTime, double dt) {
    timeSeries.Reset(startTime, endTime, dt);
  }

  void ResetTime(double startTime, double endTime, int steps) {
    timeSeries.Reset(startTime, endTime, steps);
  }

  int Step() const { return timeSeries.Step(); }

  double Time() const { return timeSeries.Time(); }

  int MaxStep() const { return timeSeries.MaxStep(); }

  double StepSize() const { return timeSeries.StepSize(); }

  double OldStepSize() const { return timeSeries.OldStepSize(); }

  bool IsFinished() const { return timeSeries.IsFinished(); }

  double LastTStep() const { return timeSeries.LastTStep(); }

  double FirstTStep() const { return timeSeries.FirstTStep(); }

  double MinStepSize() const { return timeSeries.MinStepSize(); }

  double NextTimeStep(bool iterateToNext = true) { return timeSeries.NextTimeStep(iterateToNext); }

  double NextTimeStep(double t1) { return timeSeries.NextTimeStep(t1); }

  double NextHalfTimeStep(double t1) { return timeSeries.NextHalfTStep(t1); }

  const TimeSeries &GetTimeSeries() { return timeSeries; };
};

class ILinearTimeAssemble : public ILinearAssemble, public TimeInterface {
public:
  ILinearTimeAssemble() = default;

  explicit ILinearTimeAssemble(double startTime, double endTime, double dt) :
      ILinearAssemble(), TimeInterface(startTime, endTime, dt) {}

  explicit ILinearTimeAssemble(double startTime, double endTime, int steps) :
      ILinearAssemble(), TimeInterface(startTime, endTime, steps) {}

  virtual void SetInitialValue(Vector &u) = 0;

  virtual void FinishTimeStep(const Vector &u) = 0;

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Start Time", FirstTStep(), 1),
                   PrintInfoEntry("End Time", LastTStep(), 1),
                   PrintInfoEntry("Step size", StepSize(), 1));
  }

  virtual void PrintIteration(const Vector &u) const {
    mout.PrintIteration(verbose, PrintIterEntry("n", Step(), int(log10(MaxStep())), 1),
                        PrintIterEntry("t", Time(), 7, 1), PrintIterEntry("Norm", norm(u), 13, 1));
  }

  virtual void PlotIteration(const Vector &u) const {
    // Todo
    char filename[128];
    mpp::plot("U") << u << mpp::save_plot(NumberName("U", filename, Step()));
  }
};

class INonLinearTimeAssemble : public IAssemble, public ILinearAssemble, public TimeInterface {
protected:
  std::unique_ptr<Vector> prevU;

  using IAssemble::verbose;
public:
  INonLinearTimeAssemble() = default;

  explicit INonLinearTimeAssemble(double startTime, double endTime, double dt) :
      IAssemble(), ILinearAssemble(), TimeInterface(startTime, endTime, dt) {}

  explicit INonLinearTimeAssemble(double startTime, double endTime, int steps) :
      IAssemble(), ILinearAssemble(), TimeInterface(startTime, endTime, steps) {}

  virtual void SetInitialValue(Vector &u) = 0;

  virtual void FinishTimeStep(const Vector &u) = 0;

  void Initialize(Vector &u) const override {
    THROW("Initialize(Vector &u) not implemented in INonLinearTimeAssemble")
  }

  void Initialize(Vectors &u) const override {
    THROW("Initialize(Vectors &u) not implemented in INonLinearTimeAssemble")
  }

  const char *Name() const override { return "INonLinearTimeAssemble"; }

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Start Time", FirstTStep(), 1),
                   PrintInfoEntry("End Time", LastTStep(), 1),
                   PrintInfoEntry("Step size", StepSize(), 1));
  }

  void PrintInfo(const Vector &u) const override {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Discretization", u.GetDisc().DiscName(), 2),
                   PrintInfoEntry("Level", u.GetMesh().Level().str(), 2),
                   PrintInfoEntry("Start Time", FirstTStep(), 1),
                   PrintInfoEntry("End Time", LastTStep(), 1),
                   PrintInfoEntry("Step size", StepSize(), 1));
  }

  virtual void PrintIteration(const Vector &u) const {
    mout.PrintIteration(verbose, PrintIterEntry("n", Step(), int(log10(MaxStep())), 1),
                        PrintIterEntry("t", Time(), 7, 1), PrintIterEntry("Norm", norm(u), 13, 1));
  }

  virtual void PlotIteration(const Vector &u) const {
    // Todo
    char filename[128];
    mpp::plot("U") << u << mpp::save_plot(NumberName("U", filename, Step()));
  }

  const Vector &PreviousSolution() const { return *prevU; }

  virtual void InitTimeStep(Vector &u) { prevU = std::make_unique<Vector>(u); }
};

bool isAssembleConsistent(const IAssemble &A, Vector &u);

#endif // of #ifndef _ASSEMBLE_H_
