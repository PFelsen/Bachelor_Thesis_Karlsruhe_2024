#ifndef UNIFORMDISTRIBUTION_HPP
#define UNIFORMDISTRIBUTION_HPP

#include "SampleGenerator.hpp"
#include "sprng_cpp.h"

class Random {
public:
  Random() = delete;

  static void Initialize();

  Random(const Random &random) = delete;

  Random &operator=(const Random &random) = delete;

  static double Uniform(int commSplit);

  static double Uniform(int commSplit, double leftBnd, double rightBnd);

  static RVector Uniform(int commSplit, int size, double leftBnd,
                         double rightBnd);

  static RMatrix Uniform(int commSplit, int rows, int cols, double leftBnd,
                         double rightBnd);

  static double Normal(int commSplit);

  static std::complex<double> ComplexNormal(int commSplit);

  static RVector Normal(int commSplit, int size);

  static CVector ComplexNormal(int commSplit, int size);

  static RMatrix Normal(int commSplit, int rows, int cols);

  static CMatrix ComplexNormal(int commSplit, int rows, int cols);

  static RTensor Normal(int commSplit, int fstComp,
                        int secComp, int thirdComp);

  static CTensor ComplexNormal(int commSplit, int fstComp,
                               int secComp, int thirdComp);


  static double Normal(int commSplit, double mean, double var) {
    return 0.0;
  }

private:
  static int SEED;

  static int GTYPE;

  static double firstUniform;

  static double secUniform;

  static double firstNormal;

  static double secNormal;

  static double boxMuller(int commSplit);

  static std::map<int, Sprng *> sprngMap;

  ~Random() {
    for (auto &[commSplit, sprngInstance] : sprngMap) {
      sprngInstance->free_sprng();
      delete sprngInstance;
    }
    sprngMap.clear();
  }
};

#endif  // UNIFORMDISTRIBUTION_HPP
