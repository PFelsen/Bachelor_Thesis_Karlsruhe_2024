#ifndef MLUQ_KLEXPANSIONGENERATOR_HPP
#define MLUQ_KLEXPANSIONGENERATOR_HPP
#define KLEXPANSIONGENERATOR_HPP

#include "Config.hpp"
#include "Random.hpp"
#include "SampleGenerator.hpp"
#include "SparseGridGenerator.hpp"
#include "SymmetricCovariance.hpp"
#include "TimeSeries.hpp"

class KLExpansionGenerator : public SampleGenerator<RVector> {
private:
  double beta = 1.0;
  std::string method_ = "MC";               // Methode zur Auswahl zwischen MC und Sc
  SparseGridGenerator sparseGridGenerator_; // Sparse grid Generator
  // Random randomGenerator_; //Zufallsgenerator für Monte Carlo

  std::shared_ptr<Meshes> meshes;

  int cellCount = 0;

  static int getTruncFromConfig() {
    int trunc = 5;
    Config::Get("trunc", trunc);
    return trunc;
  }

  static int getStochLevelFromConfig() {
    int init = 3;
    Config::Get("stochLevel", init);
    return init;
  }

  static bool getLocalFromConfig() {
    bool local = false;
    Config::Get("localPolynomial", local);
    return local;
  }
protected:
  void drawSample(const SampleID &id) override;
public:
  RVector sampleVec{};

  double sampleWeight = 0;

  double sigma = 1.0;

  double length = 1.0;

  double alpha = 1.0;

  struct Eigenpair {
    RVector Eigenvalues;
    RMatrix DEigenfunctions;
  };

  Eigenpair lambda;
  // KLExpansionGenerator() :SampleGenerator() {
  //     Config::Get("beta", beta);
  // }

  // Erweiterter Konstruktor zur Initialisierung der Methode und Dimension
  explicit KLExpansionGenerator() :
      SampleGenerator(),
      sparseGridGenerator_(GridDomain(KLExpansionGenerator::getTruncFromConfig()),
                           KLExpansionGenerator::getStochLevelFromConfig(), 0,
                           getLocalFromConfig()) { // TODO
    Config::Get("beta", beta);
    Config::Get("SampleMethod", method_);
    Config::Get("sigma", sigma);
    Config::Get("length", length);
    Config::Get("alpha", alpha);
  }

  RVector DrawSample(const SampleID &id, std::shared_ptr<Meshes> meshes) {
    mout.StartBlock(Name());
    vout(1) << id << endl;
    if (this->meshes == nullptr) {
      this->meshes = meshes;
      cellCount = meshes->fine().SpaceCellCountGeometry(-1);
      sampleVec.resize(cellCount);
      lambda = CovarianceSpectrum(*meshes);
    }
    drawSample(id);
    mout.EndBlock(verbose == 0);
    return sampleVec;
  }

  double GetSampleWeight() const { return sampleWeight; }

  RVector EvalSample() const override { return sampleVec; };

  Eigenpair CovarianceSpectrum(Meshes &meshes);

  Eigenpair CovarianceSpectrum(TimeSeries &ts);

  double CovKernel(double dist);

  string Name() const override { return "KLExpansionGenerator"; }
};

// class OrnsteinUhlenbeckProcess : public KLExpansionGenerator {
// protected:
//   void drawSample(const SampleID &id) override;
//
//   RVector sampleVec{};
//
// public:
//   RVector eigenvalues;
//
//   OrnsteinUhlenbeckProcess() : KLExpansionGenerator() {};
//
//   RVector EvalSample() const override {return sampleVec;};
//
//   RVector CovarianceSpectrum(int n);
//
//   void InitGenerator(std::shared_ptr<Meshes> meshes, int init);
//
//   Scalar eigenfunction(int k, double x);
//
//   string Name() const override { return "OrnsteinUhlenbeckProcess"; }
// };


#endif //MLUQ_KLEXPANSIONGENERATOR_HPP
