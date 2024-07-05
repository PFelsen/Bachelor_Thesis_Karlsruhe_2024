#include "ObservationOperator.hpp"

ObservationOperator::ObservationOperator(
    std::vector<std::pair<std::string, CoordinateType>> observationMap) {
  // This int input is supposed to be the coordinates in case of point observation
  for (const auto &pair : observationMap) {
    observations.push_back(CreateObservationUnique(pair.first, pair.second));
  }
}

RVector ObservationOperator::GetObservation(Vector &solution) {
  std::vector<double> obs;
  for (const std::unique_ptr<Observation> &measurement : observations) {
    obs.push_back(measurement->Evaluate(solution));
  }
  return RVector(obs);
}

PointObservation::PointObservation(CoordinateType coordinates) {
  if (coordinates.size() == 1) x = Point(coordinates[0]);
  if (coordinates.size() == 2) x = Point(coordinates[0], coordinates[1]);
  if (coordinates.size() == 3) x = Point(coordinates[0], coordinates[1], coordinates[2]);
}

double PointObservation::Evaluate(Vector &solution) { return *solution(x); }

Observation *CreateObservation(std::string name, CoordinateType coordinates) {
  if (name == "Point" || name == "PointObservation") { return new PointObservation(coordinates); }
  throw std::invalid_argument("Observation type " + name
                              + " not found. "
                                "Available types are: Point.");
}

std::unique_ptr<Observation> CreateObservationUnique(std::string name, CoordinateType coordinates) {
  return std::unique_ptr<Observation>(CreateObservation(name, coordinates));
}
