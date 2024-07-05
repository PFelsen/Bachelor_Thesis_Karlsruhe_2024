#include "TimeSeries.hpp"

int TimeSeries::stepSizeToSteps() const { return int(ceil((endTime - startTime) / dt)); }

double TimeSeries::stepsToStepSize() const { return double((endTime - startTime) / maxSteps); }

void TimeSeries::valuesValid() {
  if (dt > endTime) dt = endTime;
  if (dt > maxStepSize) dt = maxStepSize;
  if (dt < minStepSize) dt = minStepSize;
}

void TimeSeries::Reset(double sT, double eT, double _dt) {
  if (sT < infty) startTime = sT;
  if (eT < infty) endTime = eT;
  if (_dt < infty) {
    dt = _dt;
    maxSteps = stepSizeToSteps();
  }

  valuesValid();
  step = 0;
  t = startTime;
}

void TimeSeries::Reset(double sT, double eT, int steps) {
  if (sT < infty) startTime = sT;
  if (eT < infty) endTime = eT;
  if (steps > 0) {
    maxSteps = steps;
    dt = stepsToStepSize();
  }

  valuesValid();
  step = 0;
  t = startTime;
}

double TimeSeries::NextHaltTStep() { return NextHalfTStep(t); }

double TimeSeries::NextHalfTStep(double t1) {
  dt *= 0.5;
  if (dt < minStepSize)
    Exit("dt is smaller than minimal step size " + std::to_string(minStepSize)) t = t1 + dt;
  return t;
}

double TimeSeries::NextTStep(bool iterateToNext) { return NextTimeStep(iterateToNext); }

double TimeSeries::NextTStep(double t1) {
  while (t1 >= t - tEps)
    t += dt;
  if (t + tTol > LastTStep()) return endTime;
  return t;
}

double TimeSeries::NextTimeStep(double t1) {
  step++;
  return NextTStep(t1);
}

double TimeSeries::NextTimeStep(bool iterateToNext) {
  if (t >= endTime) return endTime;

  double nextTime = std::min(endTime, startTime + (++step) * dt);
  if (iterateToNext) { return t = nextTime; }

  step--;
  return nextTime;
}

void TimeSeries::readConfig() {
  Config::Get("startTime", startTime);
  Config::Get("t0", startTime);
  Config::Get("endTime", endTime);
  Config::Get("T", endTime);
  Config::Get("maxStepSize", maxStepSize);
  Config::Get("stepSize", dt);
  Config::Get("dt", dt);
  Config::Get("minStepSize", minStepSize);
  Config::Get("dt_min", minStepSize);
  Config::Get("maxSteps", maxSteps);
  Config::Get("steps", maxSteps);
  Config::Get("K", maxSteps);
  Config::Get("tEps", tEps);
  Config::Get("tTol", tTol);
}
