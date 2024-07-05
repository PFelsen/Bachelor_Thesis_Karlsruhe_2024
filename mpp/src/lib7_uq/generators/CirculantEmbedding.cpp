#include "CirculantEmbedding.hpp"


RVector Downsample(const RVector &fSample) {
  if (fSample.size() % 2 != 0) {
    Exit("Vector should be of size 2 * size ")
  }

  RVector cSample(fSample.size() / 2);
  for (int i = 0; i < cSample.size(); i++)
    cSample[i] = (fSample[2 * i] + fSample[2 * i + 1]) / 2;
  return cSample;
}

RMatrix Downsample(const RMatrix &fSample) {
  if (fSample.rows() % 2 != 0 || fSample.cols() % 2 != 0) {
    Exit("Matrix should be of size 2 * rows and 2 * cols")
  }

  RMatrix cSample(fSample.rows() / 2, fSample.cols() / 2);
  for (int i = 0; i < cSample.rows(); i++)
    for (int j = 0; j < cSample.cols(); j++)
      cSample[i][j] = (fSample[2 * i][2 * j] + fSample[2 * i][2 * j + 1]
                       + fSample[2 * i + 1][2 * j] + fSample[2 * i + 1][2 * j + 1]) / 4;
  return cSample;
}

RTensor Downsample(const RTensor &fSample) {//Failed test is caused by this
  if (fSample.FirstDimension() % 2 != 0 || fSample.SecondDimension() % 2 != 0 ||
      fSample.ThirdDimension() % 2 != 0) Exit(
      "Tensor dimension should be even ")
  RTensor cSample(fSample.FirstDimension() / 2, fSample.SecondDimension() / 2,
                  fSample.ThirdDimension() / 2);
  for (int i = 0; i < cSample.FirstDimension(); i++)
    for (int j = 0; j < cSample.SecondDimension(); j++)
      for (int k = 0; k < cSample.ThirdDimension(); k++)
        cSample(i, j, k) = (fSample(2 * i, 2 * j, 2 * k)
                            + fSample(2 * i, 2 * j + 1, 2 * k) + fSample(2 * i + 1, 2 * j, 2 * k)
                            + fSample(2 * i, 2 * j, 2 * k + 1) +
                            fSample(2 * i + 1, 2 * j, 2 * k + 1)
                            + fSample(2 * i, 2 * j + 1, 2 * k + 1) +
                            fSample(2 * i + 1, 2 * j + 1, 2 * k)
                            + fSample(2 * i + 1, 2 * j + 1, 2 * k + 1)) / 8;
  return cSample;
}

/*
 * 1D
 */

RVector CirculantEmbedding1D::ComputeSqrtEV() {
  RVector toepRow(cellCount);
  RVector toepCol(cellCount);
  RVector circRow = covariance.EmbeddedToeplitzMatrix(toepRow, toepCol);

  RVector eigenVal(circRow);
  FFT::RealToRealVector(circRow, eigenVal);

  for (int i = 0; i < eigenVal.size(); i++)
    if (eigenVal[i] > 0)
      eigenVal[i] = sqrt(eigenVal[i]);
    else Exit("Embedded matrix not positive definite")

  return eigenVal;
}

CVector CirculantEmbedding1D::GenerateField(const SampleID &id) {
  CVector Xf(sqrtEigenvalues);
  CVector normalDistCVector = Random::ComplexNormal(id.commSplit, sqrtEigenvalues.size());

  CVector product(sqrtEigenvalues);
  product *= normalDistCVector;

  FFT::ComplexToComplexVector(product, Xf);
  Xf /= sqrt(sqrtEigenvalues.size());
  return Xf;
}

RVector CirculantEmbedding1D::GenerateLogNormalField(const SampleID &id) {
  if (internalCounter == 0)
    fineComplexField = GenerateField(id);

  RVector logNormalField(0.0, cellCount);
  for (int i = 0; i < logNormalField.size(); i++)
    if (internalCounter == 0)
      logNormalField[i] = exp(fineComplexField[i].real() + mean);
    else
      logNormalField[i] = exp(fineComplexField[i].imag() + mean);
  internalCounter++;
  internalCounter = internalCounter % 2;
  return logNormalField;
}

RVector CirculantEmbedding1D::GenerateGaussianField(const SampleID &id) {
  if (internalCounter == 0)
    fineComplexField = GenerateField(id);

  RVector gaussianField(0.0, cellCount);
  for (int i = 0; i < gaussianField.size(); i++)
    if (internalCounter == 0)
      gaussianField[i] = fineComplexField[i].real() + mean;
    else
      gaussianField[i] = fineComplexField[i].imag() + mean;
  internalCounter++;
  internalCounter = internalCounter % 2;
  return gaussianField;
}

/*
 * 2D
 */

RMatrix CirculantEmbedding2D::ComputeSqrtEV() {
  vout(1) << "Compute Eigenvalues" << endl;
  RMatrix toepRows((int) sqrt(cellCount));
  RMatrix toepCols((int) sqrt(cellCount));
  vout(1) << "Embed Toeplitz Matrix" << endl;
  RMatrix circRows = covariance.EmbeddedToeplitzMatrix(toepRows, toepCols);

  RMatrix eigenVal(circRows);
  vout(1) << "Run FFT" << endl;
  FFT::RealToRealMatrix(circRows, eigenVal);

  for (int i = 0; i < eigenVal.rows(); i++)
    for (int j = 0; j < eigenVal.cols(); j++)
      if (eigenVal[i][j] > 0)
        eigenVal[i][j] = sqrt(eigenVal[i][j]);
      else Exit("Embedded matrix not positive definite")

  return eigenVal;
}

CMatrix CirculantEmbedding2D::GenerateField(const SampleID &id) {
  CMatrix Xf(sqrtEigenvalues);
  CMatrix normalDistCMatrix = Random::ComplexNormal(id.commSplit, sqrtEigenvalues.rows(), sqrtEigenvalues.cols());

  CMatrix product(sqrtEigenvalues);
  product *= normalDistCMatrix;

  FFT::ComplexToComplexMatrix(product, Xf);
  Xf /= sqrt(sqrtEigenvalues.size());
  return Xf;
}

RMatrix CirculantEmbedding2D::GenerateLogNormalField(const SampleID &id) {
  if (internalCounter == 0)
    fineComplexField = GenerateField(id);

  RMatrix logNormalField(0.0, (int) sqrt(cellCount));
  for (int i = 0; i < logNormalField.rows(); i++)
    for (int j = 0; j < logNormalField.cols(); j++)
      if (internalCounter == 0)
        logNormalField[i][j] = exp(fineComplexField[i][j].real() + mean);
      else
        logNormalField[i][j] = exp(fineComplexField[i][j].imag() + mean);
  internalCounter++;
  internalCounter = internalCounter % 2;
  return logNormalField;
}

RMatrix CirculantEmbedding2D::GenerateGaussianField(const SampleID &id) {
  if (internalCounter == 0)
    fineComplexField = GenerateField(id);

  RMatrix gaussianField(0.0, (int) sqrt(cellCount));
  for (int i = 0; i < gaussianField.rows(); i++)
    for (int j = 0; j < gaussianField.cols(); j++)
      if (internalCounter == 0)
        gaussianField[i][j] = fineComplexField[i][j].real() + mean;
      else
        gaussianField[i][j] = fineComplexField[i][j].imag() + mean;
  internalCounter++;
  internalCounter = internalCounter % 2;
  return gaussianField;
}

/*
 * 3D
 */

RTensor CirculantEmbedding3D::ComputeSqrtEV() {
  vout(1) << "Compute Eigenvalues" << endl;
  RTensor toepRows((int) cbrt(cellCount));
  RTensor toepCols((int) cbrt(cellCount));
  vout(1) << "Embed Toeplitz Matrix" << endl;
  RTensor circRows = covariance.EmbeddedToeplitzMatrix(toepRows, toepCols);

  RTensor eigenVal(circRows);
  vout(1) << "Run FFT" << endl;
  FFT::RealToRealTensor(circRows, eigenVal);

  for (int i = 0; i < eigenVal.FirstDimension(); i++)
    for (int j = 0; j < eigenVal.SecondDimension(); j++)
      for (int k = 0; k < eigenVal.ThirdDimension(); k++)
        if (eigenVal(i, j, k) > 0)
          eigenVal(i, j, k) = sqrt(eigenVal(i, j, k));
        else Exit("Embedded matrix not positive definite")

  return eigenVal;
}

CTensor CirculantEmbedding3D::GenerateField(const SampleID &id) {
  CTensor Xf(sqrtEigenvalues);
  CTensor normalDistCTensor =
      Random::ComplexNormal(id.commSplit, sqrtEigenvalues.FirstDimension(),
                            sqrtEigenvalues.SecondDimension(),
                            sqrtEigenvalues.ThirdDimension());

  CTensor product(sqrtEigenvalues, normalDistCTensor);

  FFT::ComplexToComplexTensor(product, Xf);
  Xf /= sqrt(sqrtEigenvalues.size());
  return Xf;
}

RTensor CirculantEmbedding3D::GenerateLogNormalField(const SampleID &id) {
  if (internalCounter == 0)
    fineComplexField = GenerateField(id);

  RTensor logNormalField(0.0, (int) cbrt(cellCount));
  for (int i = 0; i < logNormalField.FirstDimension(); i++)
    for (int j = 0; j < logNormalField.SecondDimension(); j++)
      for (int k = 0; k < logNormalField.ThirdDimension(); k++)
        if (internalCounter == 0)
          logNormalField(i, j, k) = exp(fineComplexField.real()(i, j, k) + mean);
        else
          logNormalField(i, j, k) = exp(fineComplexField.imag()(i, j, k) + mean);
  internalCounter++;
  internalCounter = internalCounter % 2;
  return logNormalField;
}

RTensor CirculantEmbedding3D::GenerateGaussianField(const SampleID &id) {
  if (internalCounter == 0)
    fineComplexField = GenerateField(id);

  RTensor gaussianField(0.0, (int) cbrt(cellCount));
  for (int i = 0; i < gaussianField.FirstDimension(); i++)
    for (int j = 0; j < gaussianField.SecondDimension(); j++)
      for (int k = 0; k < gaussianField.ThirdDimension(); k++)
        if (internalCounter == 0)
          gaussianField(i, j, k) = fineComplexField.real()(i, j, k) + mean;
        else
          gaussianField(i, j, k) = fineComplexField.imag()(i, j, k) + mean;
  internalCounter++;
  internalCounter = internalCounter % 2;
  return gaussianField;
}

Scalar CirculantEmbedding3D::EvalSample(const Point &x) const {
  int i = floor(x[0] * sample.FirstDimension());
  int j = floor(x[1] * sample.SecondDimension());
  int k = floor(x[2] * sample.ThirdDimension());

  return sample(i, j, k);
}

Scalar CirculantEmbedding3D::EvalSample(const Cell &c) const {
  return EvalSample(c());
}

void CirculantEmbedding3D::drawSample(const SampleID &id) {
  mout.StartBlock(Name());
  vout(1) << id << endl;
  if (!id.coarse) sample = GenerateLogNormalField(id);
  else sample = Downsample(sample);//sample is never initialized before?
  mout.EndBlock(verbose == 0);
}
