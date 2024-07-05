#ifndef M_CIRCULANTEMBEDDING_H
#define M_CIRCULANTEMBEDDING_H

#include <utility>

#include "SymmetricCovariance.hpp"
#include "Random.hpp"
#include "FFT.hpp"
#include "LagrangeDiscretization.hpp"

RVector Downsample(const RVector &fSample);

RMatrix Downsample(const RMatrix &fSample);

RTensor Downsample(const RTensor &fSample);


class CirculantEmbedding {

};


class CirculantEmbedding1D {
private:
  CovarianceFunction1D covariance;

  string fieldType = "LogNormal";

  std::shared_ptr<Meshes> meshes;

  CVector fineComplexField;

  RVector sqrtEigenvalues;

  int internalCounter = 0;

  double mean = 0.0;

  int cellCount = 0;

  int verbose = 0;

  RVector sample;

public:
  explicit CirculantEmbedding1D(CovarianceFunction1D covariance = CovarianceFunction1D()) :
      covariance(std::move(covariance)) {
    Config::Get("CirculantEmbeddingVerbose", verbose);
    Config::Get("StochasticField", fieldType);
    Config::Get("Mean", mean);

  }

  RVector DrawSample(const SampleID &id, std::shared_ptr<Meshes> meshes) {
    mout.StartBlock(Name());
    vout(1) << id << endl;
    if (this->meshes == nullptr) {
      this->meshes = meshes;
      cellCount = meshes->fine().SpaceCellCountGeometry(-1);
      sqrtEigenvalues = ComputeSqrtEV();
    }
    if (!id.coarse) sample = GenerateLogNormalField(id);
    else sample = Downsample(sample);
    mout.EndBlock(verbose == 0);
    return sample;
  }

  RVector ComputeSqrtEV();

  CVector GenerateField(const SampleID &id);

  RVector GenerateLogNormalField(const SampleID &id);

  RVector GenerateGaussianField(const SampleID &id);

  std::string Name() const { return "CirculantEmbedding1D"; }
};

class CirculantEmbedding2D {
private:
  CovarianceFunction2D covariance;

  string fieldType = "LogNormal";

  CMatrix fineComplexField;

  RMatrix sqrtEigenvalues;

  int internalCounter = 0;

  int cellCount = 0;

  double mean = 0.0;

  int verbose = 0;

  RMatrix sample;

  std::shared_ptr<Meshes> meshes;

  std::shared_ptr<LagrangeDiscretization> disc;

public:
  explicit CirculantEmbedding2D(CovarianceFunction2D covariance = CovarianceFunction2D())
      : covariance(std::move(covariance)) {
    Config::Get("CirculantEmbeddingVerbose", verbose);
    Config::Get("StochasticField", fieldType);
    Config::Get("Mean", mean);

  }

  RMatrix DrawSample(const SampleID &id, std::shared_ptr<Meshes> meshes) {
    mout.StartBlock(Name());
    vout(1) << id << endl;
    if (this->meshes == nullptr) {
      this->meshes = meshes;
      cellCount = meshes->fine().SpaceCellCountGeometry(-1);
      sqrtEigenvalues = ComputeSqrtEV();
    }
    if (!id.coarse) sample = GenerateLogNormalField(id);
    else sample = Downsample(sample);
    mout.EndBlock(verbose == 0);
    return sample;
  }

  Vector Draw(const SampleID &id, std::shared_ptr<Meshes> meshes) {
    mout.StartBlock(Name());
    vout(1) << id << endl;
    if (this->meshes == nullptr) {
      this->meshes = meshes;
      this->disc = std::make_shared<LagrangeDiscretization>(*meshes, 0, 1);
      cellCount = meshes->fine().SpaceCellCountGeometry(-1);
      sqrtEigenvalues = ComputeSqrtEV();
    }
    if (!id.coarse) sample = GenerateLogNormalField(id);
    else sample = Downsample(sample);

    Vector kappa(disc);
    for (cell c = kappa.cells(); c != kappa.cells_end(); ++c)
      kappa(c(), 0) = sample[floor(c()[0] * sample.rows())][floor(c()[1] * sample.cols())];

    mout.EndBlock(verbose == 0);
    return kappa;
  }




  RMatrix ComputeSqrtEV();

  CMatrix GenerateField(const SampleID &id);

  RMatrix GenerateLogNormalField(const SampleID &id);

  RMatrix GenerateGaussianField(const SampleID &id);

  std::string Name() const { return "CirculantEmbedding2D"; }
};


class CirculantEmbedding3D : public SampleGenerator<Scalar> {
private:
  CovarianceFunction3D covariance;

  string fieldType = "LogNormal";

  CTensor fineComplexField;

  RTensor sqrtEigenvalues;

  int internalCounter = 0;

  int cellCount = 0;

  double mean = 0.0;

  int verbose = 0;

  std::shared_ptr<Meshes> meshes;

  RTensor sample;

protected:
  void drawSample(const SampleID &id) override;

public:
  explicit CirculantEmbedding3D(CovarianceFunction3D covariance = CovarianceFunction3D())
      : SampleGenerator(), covariance(std::move(covariance)) {
    Config::Get("CirculantEmbeddingVerbose", verbose);
    Config::Get("StochasticField", fieldType);
    Config::Get("Mean", mean);
  }

  RTensor DrawSample(const SampleID &id, std::shared_ptr<Meshes> meshes) {
    mout.StartBlock(Name());
    vout(1) << id << endl;
    if (this->meshes == nullptr) {
      this->meshes = meshes;
      cellCount = meshes->fine().SpaceCellCountGeometry(-1);
      sqrtEigenvalues = ComputeSqrtEV();
    }
    if (!id.coarse) sample = GenerateLogNormalField(id);
    else sample = Downsample(sample);
    mout.EndBlock(verbose == 0);
    return sample;
  }

  RTensor ComputeSqrtEV();

  CTensor GenerateField(const SampleID &id);

  RTensor GenerateLogNormalField(const SampleID &id);

  RTensor GenerateGaussianField(const SampleID &id);

  Scalar EvalSample(const Cell &c) const override;

  Scalar EvalSample(const Point &x) const override;

  string Name() const override { return "CirculantEmbedding3D"; }
};


#endif //M_CIRCULANTEMBEDDING_H
