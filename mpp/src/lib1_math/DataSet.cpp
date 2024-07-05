#include "DataSet.hpp"

#include "Buffer.hpp"

DataContainer &DataContainer::operator=(const DataContainer &y) {
  data.resize(y.size());
  std::copy(y.data.begin(), y.data.end(), data.begin());
  return *this;
}

DataContainer &DataContainer::operator=(DataContainer &&y) {
  data.resize(y.size());
  std::move(y.data.begin(), y.data.end(), data.begin());
  return *this;
}

DataContainer::DataContainer(const DataContainer &other) : data(other.size()) {
  std::uninitialized_copy(other.data.begin(), other.data.end(), data.begin());
}

DataContainer::DataContainer(DataContainer &&other) noexcept : data(other.size()) {
  std::uninitialized_move(other.data.begin(), other.data.end(), data.begin());
}

// TODO: This interpolation is specific to HeartContainers and should be subclassed
DataContainer &DataContainer::operator+=(const DataContainer &y) {
  for (int i = 0; i < size(); ++i) {
    data[i] += y[i];
  }
  return *this;
}

DataContainer &DataContainer::operator-=(const DataContainer &y) {
  for (int i = 0; i < size(); ++i) {
    data[i] -= y[i];
  }
  return *this;
}

DataContainer &DataContainer::operator*=(double a) {
  for (int i = 0; i < size(); ++i) {
    data[i] *= a;
  }
  return *this;
}

DataContainer &DataContainer::operator/=(double a) {
  for (int i = 0; i < size(); ++i) {
    data[i] /= a;
  }
  return *this;
}

DataContainer::DataContainer(int n) : data(n) {
  std::uninitialized_fill(data.begin(), data.end(), 0.0);
}

double DataContainer::operator[](int i) const { return data[i]; }

int DataContainer::size() const { return data.size(); }

bool DataContainer::operator==(const DataContainer &other) const {
  if (size() != other.size()) { return false; }
  return std::equal(data.begin(), data.end(), other.data.begin());
}

void DataContainer::Set(int i, double d) { data[i] = d; }

DataContainer::DataContainer(vector<double> &&d) : data(std::move(d)) {}

Buffer &operator<<(Buffer &b, const DataContainer &d) {
  b << short(d.size());
  for (int i = 0; i < d.size(); ++i)
    b << double(d[i]);
  return b;
}

Buffer &operator>>(Buffer &b, DataContainer &d) {
  short m;
  b >> m;
  d = DataContainer(m);
  for (int i = 0; i < d.size(); ++i)
    b >> d.data[i];
  return b;
}