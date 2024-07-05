#ifndef MLUQ_OBSERVATION_OPERATOR_H
#define MLUQ_OBSERVATION_OPERATOR_H

#include "Point.hpp"
#include "ProposalGenerators.hpp"

typedef std::vector<double> CoordinateType;

class Observation {
public:
  virtual double Evaluate(Vector &solution) = 0;
};

class PointObservation : public Observation {
private:
  Point x;
public:
  PointObservation(CoordinateType coordinates);
  double Evaluate(Vector &solution);
};

class ObservationOperator {
private:
  std::vector<std::unique_ptr<Observation>> observations;
public:
  ObservationOperator(std::vector<std::pair<std::string, CoordinateType>> observationMap);
  RVector GetObservation(Vector &solution);
};

Observation *CreateObservation(std::string name, CoordinateType coordinates);
std::unique_ptr<Observation> CreateObservationUnique(std::string name, CoordinateType coordinates);

#endif // MLUQ_OBSERVATION_OPERATOR_H