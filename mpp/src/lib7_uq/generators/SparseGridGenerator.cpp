#include "SparseGridGenerator.hpp"


double GridDomain::area() const {
  if (isDefault) return pow(2.0, dim);

  double area = 1.0;
  for (int d = 0; d < dim; d++) {
    area *= abs(right()[d] - left()[d]);
  }
  return area;
}

std::vector<RVector> SparseGridGenerator::fillSampleVector() {
  std::vector<RVector> sortedPoints{};
  std::vector<double> points = grid.getPoints();
  for (auto it_start = points.begin(); it_start != points.end(); it_start += domain.dim)
    sortedPoints.emplace_back(std::vector<double>(it_start, it_start + domain.dim));
  return sortedPoints;
}

void SparseGridGenerator::drawSample(const SampleID &id) {
  if (!initialized) {
    if (useLocalPolynomialGrid) { // Different types of grids
      grid.makeLocalPolynomialGrid(domain.dim, outputs, init, 1, TasGrid::rule_localp);
    } else {
      grid.makeGlobalGrid(domain.dim, outputs, init, depth, domain.rule);
    }

    if (!domain.isDefault) {
      grid.setDomainTransform(domain.left(), domain.right());
    }

    weights = grid.getQuadratureWeights(); // Check if weights are correct -> negative values and
                                           // 3rd + 4th moment getting racked
    samples = fillSampleVector();
    mout << "Number of quadrature points:" << GetNumPoints()
         << endl; // Todo: delete and include into SLE
    initialized = true;
  }

  if (id.coarse) return;

  sample = samples[id.number];
  weight = weights[id.number];
}
