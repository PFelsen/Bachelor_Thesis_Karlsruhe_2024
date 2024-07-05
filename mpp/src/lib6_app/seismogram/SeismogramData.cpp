#include "SeismogramData.hpp"

#include "lib6_app/seismogram/base64.hpp"

std::vector<std::vector<double>> CoordAsVec(const std::vector<Point> &VEC) {
  std::vector<std::vector<double>> vecs;
  for (const auto &P: VEC) {
    vecs.push_back({P[0], P[1]});
  }
  return vecs;
}

void SeismogramData::SetMeasurement(const Point &receiver, double t, const RVector &measurement) {
  for (int i = 0; i < measurement.size(); ++i) {
    int time = ParaDat.spec.StepAtTime(t);
    int localReceiverID = (ParaDat.CoordsToLocalIndex).at(receiver);
    (*this)(localReceiverID, time, i) = measurement[i];
  }
}

double SeismogramData::GetMeasurement(const Point &receiver, double t, int index) {
    int time = ParaDat.spec.StepAtTime(t);
    int localReceiverID = (ParaDat.CoordsToLocalIndex).at(receiver);
    return (*this)(localReceiverID, time, index);
}

void SeismogramData::AddMeasurement(size_t localIndex, double t, const RVector &measurement) {
  int time = ParaDat.spec.StepAtTime(t);
  for (int i = 0; i < measurement.size(); ++i) {
    (*this)(localIndex, time, i) += measurement[i];
  }
}

void SeismogramData::AddMeasurement(const Point &receiver, double t, const RVector &measurement) {
  int time = ParaDat.spec.StepAtTime(t);
  int localIndex = (ParaDat.CoordsToLocalIndex).at(receiver);
  for (int i = 0; i < measurement.size(); ++i) {
    (*this)(localIndex, time, i) += measurement[i];
  }
}

void SeismogramData::DivideMeasurement(const Point &receiver, double t, double dividend) {
  int time = ParaDat.spec.StepAtTime(t);
  int localIndex = (ParaDat.CoordsToLocalIndex).at(receiver);
  for (int i = 0; i < ParaDat.spec.GetNrOfMeasure(); ++i) {
    (*this)(localIndex, time, i) /= dividend;
  }
}

double &
SeismogramData::operator()(size_t localReceiverID, size_t globalTimeStep, size_t component) {
  size_t localTimeStep = globalTimeStep - ParaDat.GlobalTimesOwned[localReceiverID][0];
  size_t start = ParaDat.LocalIndexToLocalRange[localReceiverID][0];
  return data[start + localTimeStep * ParaDat.spec.GetNrOfMeasure() + component];
}

const double &
SeismogramData::operator()(size_t localReceiverID, size_t globalTimeStep, size_t component) const {
  size_t localTimeStep = globalTimeStep - ParaDat.GlobalTimesOwned[localReceiverID][0];
  size_t start = ParaDat.LocalIndexToLocalRange[localReceiverID][0];
  return data[start + localTimeStep * ParaDat.spec.GetNrOfMeasure() + component];
}

void SeismogramData::WriteToFile(const std::string &filename) const {
  const ObservationSpecification &spec = ParaDat.spec;
  ExchangeBuffer exbuf;
  size_t localMaxRange=0;
  if(!ParaDat.LocalIndexToLocalRange.empty()) {
    auto maxRange =
        std::max_element(ParaDat.LocalIndexToGlobalRange.begin(),
                         ParaDat.LocalIndexToGlobalRange.end(),
                         [](auto a, auto b) { return a[1] < b[1]; });
    localMaxRange = (*maxRange)[1];
  }

  size_t TotalSize = PPM->Max(localMaxRange) + 1;

  if (!PPM->Master()) {
    if (data.size() > 0) {
      exbuf.Send(0) << ParaDat.LocalIndexToGlobalIndex << ParaDat.LocalIndexToLocalRange
                    << ParaDat.LocalIndexToGlobalRange << data.asVector();
    }
  }
  exbuf.Communicate();

  if (PPM->Master()) {
    RVector TotalData(TotalSize);
    for (int i = 0; i < ParaDat.LocalIndexToGlobalIndex.size(); ++i) {
      const auto &localRange = ParaDat.LocalIndexToLocalRange[i];
      const auto &globalRange = ParaDat.LocalIndexToGlobalRange[i];
      std::copy(data.begin() + localRange[0], data.begin() + localRange[1] ,
                TotalData.begin() + globalRange[0]);
    }
    for (int q = 0; q < PPM->Size(); ++q) {
      if (exbuf.ReceiveSize(q) == 0) continue;
      std::vector<int> LtG;
      std::vector<std::vector<int>> LtR;
      std::vector<std::vector<int>> LtS;
      std::vector<double> d;
      exbuf.Receive(q) >> LtG >> LtR >> LtS >> d;
      for (int i = 0; i < LtG.size(); ++i) {
        std::copy(d.begin() + LtR[i][0],
                  d.begin() + LtR[i][1],
                  TotalData.begin() + LtS[i][0]);
      }
    }
    nlohmann::ordered_json full;
    json ObservationSpec(nlohmann::json::value_t::object);
    ObservationSpec["Full Receiver List"] = json(CoordAsVec(spec));
    ObservationSpec["Times"] = json({spec.TimeAtStep(0),
                                     spec.GetEndTime(),
                                     spec.GetTimeDelta()});
    ObservationSpec["TimeSteps"] = json(spec.GetNumberOfTimeSteps());
    ObservationSpec["Measure Components"] = json(spec.GetNrOfMeasure());
    std::ofstream o(filename + ".json");
    json Data(nlohmann::json::value_t::object);
    const size_t receiverSize = spec.GetNrOfMeasure() * spec.GetNumberOfTimeSteps();
    for (int i = 0; i < spec.size(); ++i) {
      auto *start = TotalData.asVector().data() + i * receiverSize;
      Data[to_string(i)] = json(polfosol::b64encode(start, receiverSize * sizeof(double)));
    }
    full["Observation Specification"] = ObservationSpec;
    full["Data"] = Data;
    o << full << endl;
    o.close();
  }
  PPM->Barrier();

}

void SeismogramData::ReadFromFile(const std::string &filename) {
  if (ParaDat.LocalReceivers.empty())
    return;
  json j;
  std::ifstream ifs(filename + ".json");
  if (!ifs.good()) Exit("Failed to open file" + filename);
  nlohmann::json Seismo = nlohmann::json::parse(ifs);
  std::vector<Point> FullReceiverList;
  for (const auto &p: Seismo.at("Observation Specification").at("Full Receiver List")) {
    FullReceiverList.emplace_back(p[0], p[1]);
  }
  std::sort(FullReceiverList.begin(), FullReceiverList.end());
  for (int i = 0; i < FullReceiverList.size(); ++i) {
    if (norm((ParaDat.spec[i] - FullReceiverList[i])) > Eps) {
      THROW("Receivers of Read Seismogram are not the same as the constructed!")
    }
  }
  std::vector<double> Times{};
  for (const auto &d: Seismo.at("Observation Specification").at("Times")) {
    Times.emplace_back(d);
  }
  if (abs(Times[0] - ParaDat.spec.TimeAtStep(0)) > Eps ||
      abs(Times[1] - ParaDat.spec.GetEndTime()) > Eps ||
      abs(Times[2] - ParaDat.spec.GetTimeDelta()) > Eps) {
    THROW("Time data does not align!")
  }
  for (int i = 0; i < ParaDat.LocalReceivers.size(); ++i) {
    auto DataString = polfosol::b64decode(
        Seismo.at("Data").at(to_string(ParaDat.LocalIndexToGlobalIndex[i])));
    const auto *DataCast = reinterpret_cast<const double *>(DataString.c_str());
    std::copy(DataCast + ParaDat.GlobalTimesOwned[i][0],
              DataCast + ParaDat.GlobalTimesOwned[i][1] + 1,
              data.begin() + ParaDat.LocalIndexToLocalRange[i][0]);
  }
}

SeismogramData operator-(const SeismogramData &a, const SeismogramData &b) {
  /*if (a.GetParaDat() != b.GetParaDat()){
    THROW("SeismogramParallelizationData do not match for SeismogramData");
  }*/
  SeismogramData result = a;
  result -= b;
  return result;
}
