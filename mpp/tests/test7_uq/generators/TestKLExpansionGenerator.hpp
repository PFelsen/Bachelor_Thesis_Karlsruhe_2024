#ifndef MLUQ_TESTKLEXPANSIONGENERATOR_HPP
#define MLUQ_TESTKLEXPANSIONGENERATOR_HPP

#include "KLExpansionGenerator.hpp"
#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"

// Todo: fix/refactor

class TestKLExpansionGenerator : public Test {
protected:
  int level = 5;

  KLExpansionGenerator KLE;

  std::shared_ptr<Meshes> meshes;

  std::shared_ptr<TimeSeries> timeSeries =
      std::make_shared<TimeSeries>(TimeSeries(0, 1, abs(meshes->Level())));

  TestKLExpansionGenerator() :
      meshes(MeshesCreator("Interval").WithLevel(level).CreateShared()), KLE() {
    KLE.DrawSample(SampleID(level, 0, false), meshes);
  }
};

/*
class TestOrnsteinUhlenbeckProcess : public TestKLExpansionGenerator {
protected:
  int level = 5;

std::shared_ptr<Meshes> meshes;

OrnsteinUhlenbeckProcess OUP;
public:
TestOrnsteinUhlenbeckProcess() :
    meshes(MeshesCreator("Interval").WithLevel(level).CreateShared()), TestKLExpansionGenerator(),
    OUP() {
  OUP.InitGenerator(meshes, 200);
};
};
 */

/*
TEST_F(TestKLExpansionGenerator, Testlambdaconvg) {
  RVector eig = KLE.lambda.Eigenvalues.asVector();
  for (int i = 1; i < eig.size(); ++i) {
    double lambda1 = eig[i - 1];
    double lambda2 = eig[i];
    mout << "EV:" << i << ", " << lambda1 << ", " << lambda2 << endl;
    EXPECT_TRUE(lambda2 < lambda1);
  }
}
 */

TEST_F(TestKLExpansionGenerator, Testeigenvalues) {
  std::vector<double> lambda_analytic = {0.7388108, 0.1380038, 0.0450885, 0.0213289, 0.0122789};
  for (int i = 0; i < 5; ++i) {
    EXPECT_NEAR(KLE.lambda.Eigenvalues.asVector().at(i), lambda_analytic[i], 1e-1);
  }
}

/*
TEST_F(TestKLExpansionGenerator, TestExpectedVariance) {
  int L = 1000 * pow(2, level);
  double eigenvalue_sum = 0.0;
  double expected_var = 0.0;
  int N = timeSeries->MaxStep() + 1;

  RVector sample_mean(N);
  RVector second_mom(N);

  for (int j = 0; j < L; ++j) {
    SampleID idFine(level, j, false);
    RVector samp = KLE.DrawAndEvalSample(idFine);
    sample_mean += samp / L;
    for (int i = 0; i < N; ++i) {
      second_mom[i] += pow(samp[i], 2);
    }
  }
  for (int i = 0; i < N; ++i) {
    expected_var += (second_mom[i] - L * pow(sample_mean[i], 2)) / (L - 1);
    eigenvalue_sum += KLE.lambda.Eigenvalues.asVector()[i];
    EXPECT_NEAR(sample_mean[i], 0.0, 76.0 * (level - 1) / L);
  }
  EXPECT_NEAR(eigenvalue_sum, expected_var / N, 60.0 * (level - 1) / L);
}
 */

/*
 TEST_F(TestOrnsteinUhlenbeckProcess, Testeigenvalues) {
 std::vector<double> lambda_analytic = {0.7388108, 0.1380038, 0.0450885, 0.0213289, 0.0122789};
 for (int i = 0; i < 5; ++i) {
   EXPECT_NEAR(OUP.eigenvalues[i], lambda_analytic[i], 1e-3);
 }
}
 */

/*TEST_F(TestOrnsteinUhlenbeckProcess, Testeigenfunction){
    for (int j = 0; j < meshes->fine().CellCountGeometry(); ++j) {
        int i = 0;
        for (cell c = meshes->fine().cells(); c != meshes->fine().cells_end(); ++c) {
            double fx = OUP.eigenfunction(j,c.Center()[0]);
            double fxx = KLE.lambda.DEigenfunctions.col(j)[i];
            EXPECT_NEAR(fx, fxx, 0.1);
            i += 1;
        }
    }
}*/

// TEST_F(TestOrnsteinUhlenbeckProcess, TestExpectedVariance) {
//     int L = 1000;
//     double eigenvalue_sum = 0.0;
//     double expected_var = 0.0;
//     int N = meshes->fine().CellCountGeometry();
//
//     RVector sample_mean(200) ;
//     RVector second_mom(N);
//
//     for (int j = 0; j < L; ++j) {
//         SampleID idFine(level, j, false);
//         RVector samp = OUP.DrawAndEvalSample(idFine);
//         sample_mean += samp/L;
//         for (int i = 0; i < N; ++i) {
//             second_mom[i] += pow(samp[i],2);
//         }
//     }
//
//     for (int k = 0; k < N; ++k) {
//         expected_var += (second_mom[k] - L * pow(sample_mean[k],2))/(L-1) ;
//         eigenvalue_sum += OUP.eigenvalues[k];
//         EXPECT_NEAR(sample_mean[k], 0, 76.0*(level-1)/L);
//     }
//     EXPECT_NEAR(eigenvalue_sum, expected_var/N, 0.07);
// }

#endif // MLUQ_TESTKLEXPANSIONGENERATOR_HPP
