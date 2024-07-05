#ifndef SPARSEGRIDGENERATOR_HPP
#define SPARSEGRIDGENERATOR_HPP

#include <utility>

#if USE_TASMANIAN

#include "TasmanianSparseGrid.hpp"
#include "SampleGenerator.hpp"

// Todo:
//    * Distribute samples ob different processes

class GridDomain {
  std::pair<std::vector<double>, std::vector<double>> ab;

public:
  TasGrid::TypeOneDRule rule;

  bool isDefault = true;

  int dim;

  explicit GridDomain(int dim, TasGrid::TypeOneDRule rule = TasGrid::rule_clenshawcurtis)
      : rule(rule), dim(dim) {
    if (rule == TasGrid::rule_clenshawcurtis) {
      for (int d = 0; d < dim; d++) {
        ab.first.push_back(-1);
        ab.second.push_back(1);
      }
    } else if (rule == TasGrid::rule_localp) { //Set Domain to [0,1] for rule=localp
      for (int d = 0; d < dim; d++) {
        ab.first.push_back(0);
        ab.second.push_back(1);
      }
    } else if (rule == TasGrid::rule_gausshermite) {
      for (int d = 0; d < dim; d++) {
        ab.first.push_back(-infty);
        ab.second.push_back(infty);
      }
      // isDefault = false;
    }
  }

  explicit GridDomain(const std::pair<std::vector<double>, std::vector<double>> &ab,
                      TasGrid::TypeOneDRule rule = TasGrid::rule_clenshawcurtis) :
      ab(ab), dim(int(ab.first.size())), rule(rule) {
    if (ab.second.size() != dim) Exit("Not a domain");
    isDefault = false;
  }

  std::vector<double> left() const { return ab.first; }

  std::vector<double> right() const { return ab.second; }

  double area() const;
};

class SparseGridGenerator : public SampleGenerator<RVector> {
protected:
  TasGrid::TypeDepth depth = TasGrid::type_level;

  TasGrid::TasmanianSparseGrid grid;

  std::vector<RVector> samples{};

  std::vector<double> weights{};

  double weight = 0.0;

  RVector sample{};

  int outputs;

  int init = 3;

  GridDomain domain;

  bool useLocalPolynomialGrid = false; // Use different grid type

  std::vector<RVector> fillSampleVector();

  void drawSample(const SampleID &id) override;

public:
  explicit SparseGridGenerator(GridDomain domain, int init, int outputs = 0, bool local = false) :
      SampleGenerator(), domain(std::move(domain)), outputs(outputs), init(init),
      useLocalPolynomialGrid(local) {}

  bool initialized = false;

  double SampleWeight(const SampleID &id) const { return weights[id.number]; }

  std::vector<double> GetWeights() const { return weights; }

  std::vector<RVector> GetSamples() { return samples; }

  const GridDomain &Domain() const { return domain; }

  double SumOfWeights() const { return domain.area(); }

  RVector EvalSample() const override { return sample; }

  double EvalWeight() const { return weight; }

  int GetNumPoints() const { return grid.getNumPoints(); }

  string Name() const override { return "SparseGridGenerator"; };
};
#endif

#endif //SPARSEGRIDGENERATOR_HPP

