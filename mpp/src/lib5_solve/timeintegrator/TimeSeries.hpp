#ifndef TIMESERIES_HPP
#define TIMESERIES_HPP

#include "Config.hpp"

class TimeSeries {
private:
  double maxStepSize = infty;

  double minStepSize = 0.0;

  double startTime = 0.0;

  double endTime = 1.0;

  double tEps = 1e-15;

  double tTol = 1e-6;

  int maxSteps = 0;

  double dt = 0.0;

  double t = 0.0;

  int step = 0;

  double stepsToStepSize() const;

  int stepSizeToSteps() const;

  void valuesValid();

  void readConfig();
public:
  TimeSeries() {
    readConfig();

    if (maxSteps != 0 && dt == 0.0) this->dt = stepsToStepSize();
    else if (maxSteps == 0 && dt != 0.0) this->maxSteps = stepSizeToSteps();
    else if (maxSteps == 0 && dt == 0.0) return;
    else {
      if (maxSteps != stepSizeToSteps() && dt != stepsToStepSize()) {
        THROW("Step size dt and steps K not consistent")
      }
    }

    valuesValid();
    t = startTime;
  }

  TimeSeries(double startTime, double endTime, double dt) :
      startTime(startTime), endTime(endTime), dt(dt) {
    this->maxSteps = stepSizeToSteps();
    valuesValid();
    t = startTime;
  }

  TimeSeries(double startTime, double endTime, int steps) :
      startTime(startTime), endTime(endTime), maxSteps(steps) {
    dt = stepsToStepSize();
    valuesValid();
    t = startTime;
  }

  double NextTimeStep(bool iterateToNext = true);

  double NextHaltTStep();

  double NextTimeStep(double t);

  double NextHalfTStep(double t);

  double SecondTStep(bool m = false) { return NextTStep(FirstTStep()); }

  int MaxStep() const { return maxSteps; }

  double StepSize() const { return dt; }

  double OldStepSize() const { return dt; }

  double Time() const { return t; }

  int Step() const { return step; }

  std::string Name() { return "UniformTimeSeries"; };

  double FirstTStep() const { return startTime; };

  double LastTStep() const { return endTime; };

  double SmallStepSize() const { return minStepSize; }

  int Steps() const { return maxSteps; }

  double CurrentTStep() const { return t; };

  double NextTStep(bool iterateToNext = true);

  double NextTStep(double t1);

  double NextSmallTStep(double t1) {
    while (t1 >= t - tEps)
      t += minStepSize;
    if (t + tTol > LastTStep()) return endTime;
    return t;
  }

  double PreviousTStep(bool iterateToPrev = false) {
    if (t > startTime) {
      if (iterateToPrev) return t = std::max(startTime, t - dt);
      else return t - dt;
    } else return startTime;
  }

  bool SpecialTimeStep() { return false; }

  double MinStepSize() const { return minStepSize; }

  double MaxStepSize() const { return maxStepSize; }

  bool IsFinished() const { return (t >= endTime); }

  // Todo: should behave like Constructors
  void Reset() { Reset(infty, infty, infty); }

  void Reset(double sT, double eT, double dt);

  void Reset(double sT, double eT, int steps);
};
#endif // TIMESERIES_HPP
