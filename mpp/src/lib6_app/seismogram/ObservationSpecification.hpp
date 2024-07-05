#ifndef FWI_OBSERVATIONSPECIFICATION_HPP
#define FWI_OBSERVATIONSPECIFICATION_HPP

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include <random>
#include <cfenv>
#include "CMatrix.hpp"
#include "FFT.hpp"
#include "Parallel.hpp"
#include "Point.hpp"
#include "RVector.hpp"
#include "TimeDate.hpp"

#include "m++.hpp"

class ObservationSpecification : public std::vector<Point> {
public:
  double StartTime{};
  double EndTime{};
  double TimeDelta{};
  int NumberOfTimeSteps{};
  int SamplingFactor{1};
  int NrOfMeasureComponents{1};

  int verbose;

  void set_n_steps() {
    NumberOfTimeSteps = static_cast<int>(round((EndTime - StartTime) / TimeDelta)) + 1;
  }

public:
  ObservationSpecification() = delete;

  ObservationSpecification(const ObservationSpecification &spec)
      : vector<Point>(spec),
        StartTime(spec.StartTime),
        EndTime(spec.EndTime), TimeDelta(spec.TimeDelta){
    verbose = 0;
    Config::Get("ObservationVerbose", verbose);
    set_n_steps();
    std::sort(this->begin(), this->end());
  }

  ObservationSpecification(double start_time, double end_time, double delta_t,
                           const std::vector<Point> &receiver_positions)
      : vector<Point>(receiver_positions),
        StartTime(start_time),
        EndTime(end_time), TimeDelta(delta_t) {
    set_n_steps();
    std::sort(this->begin(), this->end());
  }

  void refine(double fac) {
    TimeDelta = TimeDelta * fac;
    set_n_steps();
  }

  [[nodiscard]] int GetSamplingFactor() const { return SamplingFactor; }

  [[nodiscard]] int GetNrOfMeasure() const { return NrOfMeasureComponents; }

  [[nodiscard]] double GetEndTime() const { return EndTime; }

  [[nodiscard]] double GetTimeDelta() const { return TimeDelta; }

  [[nodiscard]] int GetNumberOfTimeSteps() const { return NumberOfTimeSteps; }

  [[nodiscard]] double TimeAtStep(int i) const { return StartTime + i * TimeDelta; }

  [[nodiscard]] int StepAtTime(double t) const {
    return static_cast<int>(std::round((t - StartTime) / TimeDelta));
  }

  // -- IO stuff --
  void printInfo(const char *offset = "") const {
    vout(1) << offset << "FWIObservationSpecification:" << endl;
    vout(1) << offset << "  times: [" << StartTime << ", " << EndTime
         << "],  delta_t: " << TimeDelta << endl;
    vout(1) << offset << "  receiver count: " << size() << endl;
    vout(1) << offset << "  receiver positions:" << endl;
    for (const auto &r: *this)
      vout(1) << offset << "    " << r << endl;
  }

};

inline void ReadWaypoints(const string &Wpname, int &WpCount, vector<Point> &WpPoints,
                          std::vector<int> &nrPointsBetween) {
  int tmp = WpCount;
  Config::Get(Wpname + "Waypoints", WpCount);
  if (WpCount == 0) return;
  if (tmp == WpCount) {
    if (nrPointsBetween.empty()) nrPointsBetween = std::vector<int>(WpCount - 1, 0);
    if (WpPoints.empty()) WpPoints = std::vector<Point>(WpCount, Infty);
  } else {
    size_t sizeold = WpPoints.size();
    nrPointsBetween.resize(WpCount - 1);
    WpPoints.resize(WpCount);
    if (sizeold < WpCount) {
      for (size_t i = sizeold; i < WpCount; ++i) {
        WpPoints.at(i) = Infty;
      }
      for (size_t i = sizeold; i < WpCount - 1; ++i) {
        nrPointsBetween.at(i) = 0;
      }
    }
  }
  for (int i = 0; i < WpCount; ++i) {
    // mout << (Wpname + "Count" + std::to_string(i)).c_str() << endl;
    if (i != WpCount - 1) Config::Get(Wpname + "Count" + std::to_string(i), nrPointsBetween[i]);
    if (i < WpCount) {
      Config::Get(Wpname + "Waypoint" + std::to_string(i), WpPoints[i]);
      if (WpPoints[i] == Infty) THROW("Failed to specify ReceiverWaypoint " + std::to_string(i));
    }
  }
  if (WpCount < 2) Exit("At least 2 waypoints needed!" + Wpname);
}

inline void printResult(double t0, double dt, double T, std::vector<Point> &receiver_positions) {
  int verbose = 0;
  Config::Get("ObservationVerbose", verbose, true);
  vout(1) << "ReceiverStart: " << t0 << " ";
  vout(1) << "ReceiverEnd: " << T << " ";
  vout(1) << "ReceiverDt: " << dt << endl;
  int i = 0;
  for (auto p = receiver_positions.begin(); p != receiver_positions.end(); ++p, ++i) {
    vout(1) << "  " << i << ": " << *p << " ";
    if (i < 10) vout(1) << " ";
    if (((i + 1) % p->SpaceDim()) == 0) vout(1) << endl;
  }
  vout(1) << endl;
}

inline std::vector<Point> PointsOnLine(vector<Point> &waypoints,
                                       vector<int> &nrOfPointsOnLine) {
  size_t waypoint_count = waypoints.size();
  std::vector<Point> PointsOnLine;
  Point start = waypoints[0];
  for (int i = 1; i < waypoint_count; ++i) {
    Point end = waypoints[i];
    int count = nrOfPointsOnLine[i - 1];
    if (count == 1 && norm(start - end) < Eps) {
      PointsOnLine.push_back(end);
      continue;
    }
    if (count == 1 && norm(start - end) > Eps) {
      THROW("Only one point on line but start and end are not the same!")
    }
    for (int j = 0; j < count; ++j) {
      double lambda = double(j) / (count - 1);
      PointsOnLine.push_back((1.0 - lambda) * start + lambda * end);
    }
    start = end;
  }
  return PointsOnLine;
}

inline std::vector<Point> GetReceiverPositionsOnLine(
    vector<Point> &waypoints, vector<int> &receiverCounts) {
  auto waypoint_count = waypoints.size();
  if (waypoint_count - 1 != receiverCounts.size()) Exit("Specify one count fewer than waypoints");
  if (waypoint_count < 2) Exit("At least 2 waypoints needed!");
  std::vector<Point> receiver_positions = PointsOnLine(waypoints, receiverCounts);
  return receiver_positions;
}


class ObservationSpecificationBuilder {
  double t0 = 0.0;
  double dt = -1;
  double T = -1;
  int WaypointCount = 0;
  std::vector<int> receiverCounts{};
  std::vector<Point> waypoints{};
  std::string fn;
  std::vector<Point> receivers{};
  string readStyle = "Classic";

  int verbose = 0;

  void readConfig() {
    Config::Get("ObservationVerbose", verbose, true);
    Config::Get("Specfile", fn);
    Config::Get("ReceiverReadStyle", readStyle);
    if (!fn.empty() && readStyle != "Classic") {
      THROW("More than one reading option specified!");
    }
    if (readStyle == "Points") {
      size_t nr = 0;
      Config::Get("NrReceivers", nr);
      for (int i = 0; i < nr; ++i) {
        Point tmp;
        Config::Get("Receiver" + to_string(i), tmp);
        receivers.push_back(tmp);
      }
      return;
    } else if (readStyle != "Classic") {
      THROW("Unknown ReadStyle!")
    }
    if (!fn.empty()) return;
    Config::Get("t0", t0);
    Config::Get("dt", dt);
    Config::Get("T", T);
    ReadWaypoints("Receiver", WaypointCount, waypoints, receiverCounts);
  }

public:
  ObservationSpecificationBuilder() {
    readConfig();
  }

  ObservationSpecificationBuilder& WithT(double T_) {
    T = T_;
    return *this;
  }

  ObservationSpecificationBuilder& Withdt(double dt_) {
    dt = dt_;
    return *this;
  }

  ObservationSpecificationBuilder& Witht0(double t0_) {
    t0 = t0_;
    return *this;
  }

  ObservationSpecificationBuilder& WithWaypoints(const std::vector<Point> &waypoints_) {
    waypoints = waypoints_;
    return *this;
  }

  ObservationSpecificationBuilder& WithReceiverCounts(const std::vector<int> &receiverCounts_) {
    receiverCounts = receiverCounts_;
    return *this;
  }
  ObservationSpecificationBuilder &WithReceivers(
      const std::vector<Point> &receivers_) {
    receivers = receivers_;
    return *this;
  }

  ObservationSpecification *Create() {
    if (fn.empty()) {
      if (receivers.empty() && !waypoints.empty() && !receiverCounts.empty()) {
        receivers = GetReceiverPositionsOnLine(waypoints, receiverCounts);
      } else if (receivers.empty() && waypoints.empty() &&
                 receiverCounts.empty()) {
        THROW("Specify either Specfile, waypoints or explicit receivers in "
              "creation");
      }
      return new ObservationSpecification(t0, T, dt, receivers);
    }
    return CreateFromFile();
  }

  std::shared_ptr<ObservationSpecification> CreateShared() {
    return std::shared_ptr<ObservationSpecification>(Create());
  }

  std::unique_ptr<ObservationSpecification> CreateUnique() {
    return std::unique_ptr<ObservationSpecification>(Create());
  }

  ObservationSpecification *CreateFromFile() {
    json j;
    ExchangeBuffer exbuf;
    std::ifstream ifs(fn + ".json");
    if (!ifs.good()) Exit("Failed to open file" + fn + ".json");
    nlohmann::json Seismo = nlohmann::json::parse(ifs);
    std::vector<Point> FullReceiverList;
    for (const auto &p: Seismo.at("Observation Specification").at("Full Receiver List")) {
      FullReceiverList.emplace_back(p[0], p[1]);
    }
    std::sort(FullReceiverList.begin(), FullReceiverList.end());
    std::vector<double> Times{};
    for (const auto &d: Seismo.at("Observation Specification").at("Times")) {
      Times.emplace_back(d);
    }
    return new ObservationSpecification(Times[0], Times[1], Times[2], FullReceiverList);
  }
};

#endif  // FWI_OBSERVATIONSPECIFICATION_HPP
