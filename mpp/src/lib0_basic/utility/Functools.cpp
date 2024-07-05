#include "Functools.hpp"
#include "Assertion.hpp"

std::vector<double> ft::linspace(double init, double end, unsigned int size) {
  if (size < 2) THROW("Not Defined for size < 2");
  std::vector<double> lin;
  double inc = (end - init) / (size - 1);
  std::generate_n(std::back_inserter(lin), size,
                  [n = 0.0, init, inc]() mutable { return init + (n++) * inc; });
  return lin;
}

std::pair<double, double> ft::linearFit(std::vector<double> &x, std::vector<double> &y) {
  long xSize = x.size();
  double xSum = 0, ySum = 0, xxSum = 0, xySum = 0, slope, intercept;
  for (long i = 0; i < xSize; i++) {
    xSum += x[i];
    ySum += y[i];
    xxSum += x[i] * x[i];
    xySum += x[i] * y[i];
  }
  slope = ((double)xSize * xySum - xSum * ySum) / ((double)xSize * xxSum - xSum * xSum);
  intercept = (ySum - slope * xSum) / ((double)xSize);
  std::pair<double, double> res;
  res.first = slope;
  res.second = intercept;
  return res;
}