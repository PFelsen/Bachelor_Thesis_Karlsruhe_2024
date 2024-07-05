#include "Algebra.hpp"
#include "CirculantEmbedding.hpp"
#include "LagrangeDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "Plotting.hpp"

#include <filesystem>
#include "TestEnvironment.hpp"

namespace fs = std::filesystem;

class TestCirculantEmbedding1D : public Test {
protected:
  std::shared_ptr<Meshes> meshes;

  CovarianceFunction1D covariance;

  CirculantEmbedding1D circEmb;

  int level = 7;

  int samples = 2;

  int commSplit;

  explicit TestCirculantEmbedding1D(
      int commSplit = 0, CovarianceFunction1D covariance = CovarianceFunction1D(1.0, 0.15, 1.8)) :
      covariance(covariance), circEmb(covariance), commSplit(commSplit) {

    meshes = MeshesCreator("Interval").WithPLevel(level - 1).WithLevel(level).CreateShared();
  }

  void TearDown() override { PPM->Barrier(0); }
};

class TestCirculantEmbedding1DWithSplit : public TestCirculantEmbedding1D {
public:
  TestCirculantEmbedding1DWithSplit() : TestCirculantEmbedding1D(1){};
};

class TestCirculantEmbedding1DWithDoubleSplit : public TestCirculantEmbedding1D {
public:
  TestCirculantEmbedding1DWithDoubleSplit() : TestCirculantEmbedding1D(2){};
};

class TestCirculantEmbedding2D : public Test {
protected:
  std::shared_ptr<Meshes> meshes;

  CovarianceFunction2D covariance;

  CirculantEmbedding2D circEmb;

  int level = 5;

  int samples = 2;

  int commSplit;

  explicit TestCirculantEmbedding2D(
      int commSplit = 0,
      CovarianceFunction2D covariance = CovarianceFunction2D(1.0, {0.15, 0.15}, 1.8, 2)) :
      covariance(covariance), circEmb(covariance), commSplit(commSplit) {

    meshes = MeshesCreator("Square").WithPLevel(level - 1).WithLevel(level).CreateShared();
  }

  void TearDown() override { PPM->Barrier(0); }
};

class TestCirculantEmbedding2DWithSplit : public TestCirculantEmbedding2D {
public:
  TestCirculantEmbedding2DWithSplit() : TestCirculantEmbedding2D(1){};
};

class TestCirculantEmbedding2DWithDoubleSplit : public TestCirculantEmbedding2D {
public:
  TestCirculantEmbedding2DWithDoubleSplit() : TestCirculantEmbedding2D(2){};
};

class TestCirculantEmbedding3D : public Test {
protected:
  std::shared_ptr<Meshes> meshes;

  CovarianceFunction3D covariance;

  CirculantEmbedding3D circEmb;

  int level = 6;

  int samples = 2;

  int commSplit;

  explicit TestCirculantEmbedding3D(
      int commSplit = 0,
      CovarianceFunction3D covariance = CovarianceFunction3D(1.0, {0.15, 0.15, 0.15}, 1.8, 2)) :
      covariance(covariance), circEmb(covariance), commSplit(commSplit) {

    meshes = MeshesCreator("Hexahedron").WithPLevel(level - 1).WithLevel(level).CreateShared();
  }

  void TearDown() override { PPM->Barrier(0); }
};

class TestCirculantEmbedding3DWithSplit : public TestCirculantEmbedding3D {
public:
  TestCirculantEmbedding3DWithSplit() : TestCirculantEmbedding3D(1){};
};

class TestCirculantEmbedding3DWithDoubleSplit : public TestCirculantEmbedding3D {
public:
  TestCirculantEmbedding3DWithDoubleSplit() : TestCirculantEmbedding3D(2){};
};

#define CIRCULANTEMBEDDING_DRAWSAMPLE(TestClass)                                                   \
                                                                                                   \
  TEST_F(TestClass, TestDrawSample) {                                                              \
    for (int sample = 0; sample < samples; sample++) {                                             \
                                                                                                   \
      SampleID idFine(level, sample, false, "Kappa");                                              \
      idFine = idFine.WithCommSplit(commSplit);                                                    \
      auto y = circEmb.DrawSample(idFine, meshes);                                                 \
      SampleID idCoarse(level, sample, true, "Kappa");                                             \
      idCoarse = idFine.WithCommSplit(commSplit);                                                  \
      circEmb.DrawSample(idCoarse, meshes);                                                        \
    }                                                                                              \
  }

// Todo: Think about TestCase that tries to estimate mean and variance.
//       Consider plots for 3D

class Plot2DSmoothing10Lambda005P1 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing10Lambda005P1() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.05, 0.05}, 1.0, 1)) {}
};

class Plot2DSmoothing10Lambda010P1 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing10Lambda010P1() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.10, 0.05}, 1.0, 1)) {}
};

class Plot2DSmoothing10Lambda015P1 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing10Lambda015P1() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.15, 0.05}, 1.0, 1)) {}
};

class Plot2DSmoothing10Lambda005P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing10Lambda005P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.05, 0.05}, 1.0, 2)) {}
};

class Plot2DSmoothing10Lambda010P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing10Lambda010P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.10, 0.05}, 1.0, 2)) {}
};

class Plot2DSmoothing10Lambda015P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing10Lambda015P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.15, 0.05}, 1.0, 2)) {}
};

class Plot2DSmoothing14Lambda005P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing14Lambda005P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.05, 0.05}, 1.4, 2)) {}
};

class Plot2DSmoothing14Lambda010P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing14Lambda010P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.10, 0.05}, 1.4, 2)) {}
};

class Plot2DSmoothing14Lambda015P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing14Lambda015P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.15, 0.05}, 1.4, 2)) {}
};

class Plot2DSmoothing18Lambda005P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing18Lambda005P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.05, 0.05}, 1.8, 2)) {}
};

class Plot2DSmoothing18Lambda010P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing18Lambda010P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.10, 0.05}, 1.8, 2)) {}
};

class Plot2DSmoothing18Lambda015P2 : public TestCirculantEmbedding2D {
public:
  Plot2DSmoothing18Lambda015P2() :
      TestCirculantEmbedding2D(0, CovarianceFunction2D(1.0, {0.15, 0.05}, 1.8, 2)) {}
};

/*
 * 3D
 */

class Plot3DSmoothing10Lambda015P2 : public TestCirculantEmbedding3D {
public:
  Plot3DSmoothing10Lambda015P2() :
      TestCirculantEmbedding3D(0, CovarianceFunction3D(1.0, {0.15, 0.15, 0.15}, 1.0, 2)) {}
};

class Plot3DSmoothing10Lambda010P2 : public TestCirculantEmbedding3D {
public:
  Plot3DSmoothing10Lambda010P2() :
      TestCirculantEmbedding3D(0, CovarianceFunction3D(1.0, {0.15, 0.05, 0.10}, 1.0, 2)) {}
};

class Plot3DSmoothing10Lambda005P2 : public TestCirculantEmbedding3D {
public:
  Plot3DSmoothing10Lambda005P2() :
      TestCirculantEmbedding3D(0, CovarianceFunction3D(1.0, {0.05, 0.05, 0.05}, 1.0, 2)) {}
};

class Plot3DSmoothing18Lambda015P2 : public TestCirculantEmbedding3D {
public:
  Plot3DSmoothing18Lambda015P2() :
      TestCirculantEmbedding3D(0, CovarianceFunction3D(1.0, {0.15, 0.15, 0.15}, 1.8, 2)) {}
};

class Plot3DSmoothing10Lambda015P1 : public TestCirculantEmbedding3D {
public:
  Plot3DSmoothing10Lambda015P1() :
      TestCirculantEmbedding3D(0, CovarianceFunction3D(1.0, {0.15, 0.15, 0.15}, 1.0, 1)) {}
};

#define CIRCULANTEMBEDDING_PLOTTING(TestClass)                                                     \
                                                                                                   \
  TEST_F(TestClass, TestPlotting) {                                                                \
    /* TODO: This always uses commsplit 0 */                                                       \
    SampleID idFine(level, 0, false, covariance.ToString());                                       \
    idFine = idFine.WithCommSplit(commSplit);                                                      \
    Vector kappa = circEmb.Draw(idFine, meshes);                                                   \
    mpp::plot(covariance.ToString()) << kappa << mpp::endp;                                        \
    std::string filePVTU = Config::GetPlotPath() + covariance.ToString() + ".pvtu";                \
    std::string fileVTU = Config::GetPlotPath() + covariance.ToString() + ".vtu";                  \
    PPM->Barrier(0);                                                                               \
    EXPECT_TRUE(fs::exists(filePVTU) || fs::exists(fileVTU));                                                                 \
  }


CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding1D)

CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding1DWithSplit)

CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding1DWithDoubleSplit)

#if (SpaceDimension >= 2)
CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding2D)

CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding2DWithSplit)

CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding2DWithDoubleSplit)
#endif

#if (SpaceDimension == 3)
// CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding3D)
//
// CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding3DWithSplit)
//
// CIRCULANTEMBEDDING_DRAWSAMPLE(TestCirculantEmbedding3DWithDoubleSplit)
#endif

/*
 * PlotTests 2D - 3D
 */

#if (SpaceDimension >= 2)
CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing10Lambda005P1)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing10Lambda010P1)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing10Lambda015P1)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing10Lambda005P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing10Lambda010P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing10Lambda015P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing14Lambda005P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing14Lambda010P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing14Lambda015P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing18Lambda005P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing18Lambda010P2)

CIRCULANTEMBEDDING_PLOTTING(Plot2DSmoothing18Lambda015P2)
#endif

#if (SpaceDimension == 3)
// CIRCULANTEMBEDDING_PLOTTING(Plot3DSmoothing10Lambda015P2)
//
// CIRCULANTEMBEDDING_PLOTTING(Plot3DSmoothing10Lambda005P2)
//
// CIRCULANTEMBEDDING_PLOTTING(Plot3DSmoothing18Lambda015P2)
//
// CIRCULANTEMBEDDING_PLOTTING(Plot3DSmoothing10Lambda015P1)
#endif

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("GeneratorVerbose", 0)
                     .WithScreenLogging()
                     .WithRandomInitialized()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}