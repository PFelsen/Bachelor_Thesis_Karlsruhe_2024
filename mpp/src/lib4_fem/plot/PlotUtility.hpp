#ifndef PLOTUTILITY_HPP
#define PLOTUTILITY_HPP

#include <vector>

enum PlotStatus { CLEAR, FLUSH, SAVE };

class Buffer;

class PointData {
  std::vector<double> data;
public:
  explicit PointData(int dim = 1) : data(dim) {}

  explicit PointData(const std::vector<double> &d) : data(d) {}

  explicit PointData(std::vector<double> &&d) : data(std::move(d)) {}

  void AddData(int index, double d) { data[index] = d; }

  double GetData(int i) const { return data.at(i); }

  // TODO: Research vector.back()
  double PlotData() { return data.back(); }

  friend Buffer &operator<<(Buffer &buffer, const PointData &pointData);
  friend Buffer &operator>>(Buffer &buffer, PointData &pointData);
};

#endif // PLOTUTILITY_HPP
