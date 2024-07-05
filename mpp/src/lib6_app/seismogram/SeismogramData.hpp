#ifndef FWI_SEISMOGRAMDATA_HPP
#define FWI_SEISMOGRAMDATA_HPP


#include "ObservationSpecification.hpp"
#include "Parallel.hpp"
#include "Point.hpp"
#include "RVector.hpp"
#include "TimeDate.hpp"

#include "m++.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <utility>

class SeismogramParallelizationData {
public:
  std::unordered_map<Point, int> CoordsToLocalIndex;
  std::vector<Point> LocalReceivers;
  vector<int> LocalIndexToGlobalIndex{};
  std::vector<std::vector<int>> LocalIndexToLocalRange;
  std::vector<std::vector<int>> LocalIndexToGlobalRange;
  std::vector<std::vector<int>> GlobalTimesOwned;
  const ObservationSpecification &spec;

  explicit SeismogramParallelizationData(const ObservationSpecification &spec_) : spec(spec_) {
  }

};

class SeismogramData {
  const SeismogramParallelizationData ParaDat;
  RVector data{};
public:
  SeismogramData() = delete;

  SeismogramData(const SeismogramData &s) : ParaDat(s.GetParaDat()) {
    data = s.data;
  }

  explicit SeismogramData(SeismogramParallelizationData ParaDat_) : ParaDat(std::move(ParaDat_)) {
    int size = 0;
    for (const auto &LtR: ParaDat.LocalIndexToLocalRange) {
      size = max(size, LtR[1]);
    }
    data = RVector(0.0, size + 1);
  }

  void SetMeasurement(const Point &receiver, double t, const RVector &measurement);

  double GetMeasurement(const Point &receiver, double t, int index);

  void AddMeasurement(const Point &receiver, double t, const RVector &measurement);

  void AddMeasurement(size_t localIndex, double t, const RVector &measurement);

  void DivideMeasurement(const Point &receiver, double t, double dividend);

  double &operator()(size_t localReceiverID, size_t globalTimeStep, size_t component);

  const double &operator()(size_t localReceiverID, size_t globalTimeStep, size_t component) const;

  SeismogramData &operator=(double val) {
    data = val;
    return *this;
  }

  SeismogramData &operator=(const SeismogramData &dat) {
    data = dat.data;
    return *this;
  }

  SeismogramData &operator=(const SeismogramData &&dat) {
    data = dat.data;
    return *this;
  }

  void operator*=(double val) {
    data *= val;
  }

  void operator-=(const SeismogramData &s) {
    data -= s.data;
  }

  void operator+=(const SeismogramData &s) {
    data += s.data;
  }

  const RVector &GetData() const { return data; }

  const SeismogramParallelizationData &GetParaDat() const {
    return ParaDat;
  }

  std::vector<int> GetLocToRange(size_t i) const {
    return ParaDat.LocalIndexToLocalRange[i];
  }

  size_t ReceiverSize() const { return ParaDat.LocalReceivers.size(); }

  auto GetTimeseries(const Point &receiver) const {
    const auto &localIndex = ParaDat.CoordsToLocalIndex.at(receiver);
    const auto &localRange = ParaDat.LocalIndexToLocalRange[localIndex];
    auto start = data.asVector().begin();
    return std::span{start + localRange[0], start + localRange[1] + 1};
  }

  auto GetTimeseries(const Point &receiver) {
    const auto &localIndex = ParaDat.CoordsToLocalIndex.at(receiver);
    const auto &localRange = ParaDat.LocalIndexToLocalRange[localIndex];
    auto start = data.asVector().begin();
    return std::span{start + localRange[0], start + localRange[1] + 1};
  }

  void WriteToFile(const std::string &filename) const;

  void ReadFromFile(const std::string &filename);
};

SeismogramData operator-(const SeismogramData &a, const SeismogramData &b);


#endif  // FWI_SEISMOGRAMDATA_HPP
