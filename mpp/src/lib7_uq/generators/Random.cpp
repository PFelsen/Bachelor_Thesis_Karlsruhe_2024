#include <cmath>

#include "Random.hpp"

int Random::GTYPE = 1;

int Random::SEED = 985456376;

double Random::firstUniform = 0.0;

double Random::secUniform = 0.0;

double Random::firstNormal = 0.0;

double Random::secNormal = 0.0;

std::map<int, Sprng *> Random::sprngMap = {};

void Random::Initialize() {
  if (sprngMap.empty()) {
    for (int commSplit = 0; commSplit <= PPM->MaxCommSplit(); commSplit++) {
      sprngMap[commSplit] = SelectType(GTYPE);
      sprngMap[commSplit]->init_sprng(
          PPM->Color(commSplit), (int) pow(2, commSplit), SEED + commSplit, SPRNG_DEFAULT
      );
    }
  }
}

double Random::Uniform(int commSplit) {
  return sprngMap[commSplit]->sprng();
}

double Random::Uniform(int commSplit, double leftBnd, double rightBnd) {
  return leftBnd + Random::Uniform(commSplit) * (rightBnd - leftBnd);
}

RVector Random::Uniform(int commSplit, int size, double leftBnd,
                        double rightBnd) {
  RVector results(size);
  for (int i = 0; i < size; i++)
    results[i] = Random::Uniform(commSplit, leftBnd, rightBnd);
  return results;
}

RMatrix Random::Uniform(int commSplit, int rows, int cols, double leftBnd,
                        double rightBnd) {
  RMatrix results(rows, cols);
  for (int i = 0; i < rows; i++) {
    results.InsertRow(Random::Uniform(commSplit, cols, leftBnd, rightBnd), i);
  }
  return results;
}

double Random::boxMuller(int commSplit) {
  if (firstUniform == 0.0 && secUniform == 0.0) {
    firstUniform = Random::Uniform(commSplit);
    secUniform = Random::Uniform(commSplit);
    firstNormal = sqrt(-2.0 * log(firstUniform)) * cos(2.0 * M_PI * secUniform);
    secNormal = sqrt(-2.0 * log(firstUniform)) * sin(2.0 * M_PI * secUniform);
    return firstNormal;
  } else {
    firstUniform = 0.0;
    secUniform = 0.0;
    return secNormal;
  }
}

double Random::Normal(int commSplit) {
  return boxMuller(commSplit);
}

std::complex<double> Random::ComplexNormal(int commSplit) {
  return {Random::Normal(commSplit), Random::Normal(commSplit)};
}

RVector Random::Normal(int commSplit, int size) {
  RVector results(size);
  for (int i = 0; i < size; i++)
    results[i] = Random::Normal(commSplit);
  return results;
}

CVector Random::ComplexNormal(int commSplit, int size) {
  CVector results(size);
  for (int i = 0; i < size; i++)
    results[i] = Random::ComplexNormal(commSplit);
  return results;
}

RMatrix Random::Normal(int commSplit, int rows, int cols) {
  RMatrix results(rows, cols);
  for (int i = 0; i < rows; i++) {
    results.InsertRow(Random::Normal(commSplit, cols), i);
  }
  return results;
}

CMatrix Random::ComplexNormal(int commSplit, int rows, int cols) {
  CMatrix results(rows, cols);
  for (int i = 0; i < rows; i++) {
    results.InsertRow(Random::ComplexNormal(commSplit, cols), i);
  }
  return results;
}

RTensor Random::Normal(int commSplit, int fstComp, int secComp, int thirdComp) {
  RTensor results(fstComp, secComp, thirdComp);
  for (int i = 0; i < thirdComp; i++) {
    results.InsertMatrixThirdDimension(
        Random::Normal(commSplit, fstComp, secComp), i
    );
  }
  return results;
}

CTensor Random::ComplexNormal(int commSplit, int fstComp, int secComp, int thirdComp) {
  CTensor results(fstComp, secComp, thirdComp);
  for (int i = 0; i < thirdComp; i++) {
    results.InsertMatrixThirdDimension(
        Random::ComplexNormal(commSplit, fstComp, secComp), i
    );
  }
  return results;
}
