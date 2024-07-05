#ifndef _DATASET_H_
#define _DATASET_H_

#include <utility>

#include "Tensor.hpp"

/**
 * \class DataContainer
 * \brief Stores data from an appropriate .geo file.
 *
 * DataContainers are used when geometries with additional vertex data are loaded.
 * Arithmetic operators are needed for mesh refinements, as the data on the newly created
 * vertices has to be calculated from the coarse mesh.
 */
class DataContainer {
public:
  std::vector<double> data{};

  explicit DataContainer(int n = 1);

  explicit DataContainer(std::vector<double> &&d);

  DataContainer(const DataContainer &other);

  DataContainer(DataContainer &&other) noexcept;

  DataContainer &operator=(const DataContainer &y);

  DataContainer &operator=(DataContainer &&y);

  double operator[](int i) const;

  /// Sets the value at index i to the value d.
  void Set(int i, double d);

  /// Returns size of this DataContainer.
  int size() const;

  DataContainer &operator+=(const DataContainer &y);

  DataContainer &operator-=(const DataContainer &y);

  DataContainer &operator*=(double a);

  DataContainer &operator/=(double a);

  bool operator==(const DataContainer &other) const;

  auto begin() { return data.begin(); }

  auto end() { return data.end(); }
};

inline DataContainer operator+(const DataContainer &x, const DataContainer &y) {
  DataContainer z = x;
  return z += y;
}

inline DataContainer operator-(const DataContainer &x, const DataContainer &y) {
  DataContainer z = x;
  return z -= y;
}

inline DataContainer operator*(double b, const DataContainer &x) {
  DataContainer z = x;
  return z *= b;
}

inline DataContainer operator*(const DataContainer &x, double b) {
  DataContainer z = x;
  return z *= b;
}

inline DataContainer operator*(int b, const DataContainer &x) {
  DataContainer z = x;
  return z *= double(b);
}

inline DataContainer operator/(const DataContainer &x, double b) {
  DataContainer z = x;
  return z /= b;
}

inline std::ostream &operator<<(std::ostream &os, const DataContainer &z) {
  os << "DataContainer (" << z.size() << "):";
  for (int i = 0; i < z.size(); ++i)
    os << " - " << z[i];
  return os;
}

class Buffer;

Buffer &operator<<(Buffer &b, const DataContainer &d);

Buffer &operator>>(Buffer &b, DataContainer &d);


#endif
