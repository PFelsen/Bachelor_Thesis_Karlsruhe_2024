#include "TestEnvironment.hpp"
#include "Point.hpp" // Replace with the correct header file for your Point class
#include "KLExpansion.hpp" // Replace with the correct header for your SimpleKLExpansion class
#include "RVector.hpp" // Replace with the correct header for your RVector class

TEST(PointTest, InitializationAndAccess) {
    double coordinate1 = 0.5;
    double coordinate2 = 0.7;
    Point x(coordinate1, coordinate2);

    EXPECT_DOUBLE_EQ(x[0], coordinate1);
    EXPECT_DOUBLE_EQ(x[1], coordinate2);
}

TEST(SimpleKLExpansionTest, EvaluateSample) {
    RVector randomInput(1); 
    randomInput[0] = 1;

    SimpleKLExpansion klExpansion;
    klExpansion.GetParameterVector() = randomInput;
    Point testPoint(0.5, 0.5);

    EXPECT_DOUBLE_EQ(klExpansion.EvalSample(testPoint), 1);
}

TEST(TwoDimKLExpansionTest, EvaluateSample) {
    RVector randomInput(0, 4);
    randomInput[0] = 1;      

    TwoDimKLExpansion klExpansion;
    klExpansion.GetParameterVector() = randomInput;
    Point testPoint(0.5, 0);

    double expected = exp(0.15 + std::sqrt(0.03));
//TODO: FIX ME
//    EXPECT_NEAR(klExpansion.EvalSample(testPoint), expected, 1e-6);
}

TEST(ComplexTwoDimKLExpansionTest, EvaluateSample) {
    RVector randomInput(0, 64);
    randomInput[0] = 1;      

    ComplexTwoDimKLExpansion klExpansion;
    klExpansion.GetParameterVector() = randomInput;
    Point testPoint(0, 0);

    double expected = exp(0.15 + std::sqrt(0.1));

    EXPECT_NEAR(klExpansion.EvalSample(testPoint), expected, 1e-6);
}

TEST(CreateRandomFieldUniqueTest, CorrectClassCreation) {
    RVector randomInput;

    auto simpleKLExpansion = CreateRandomFieldUnique("Simple2DKL");
    EXPECT_EQ(simpleKLExpansion->Name(), "Simple2DKLExpansion");

    auto twoDimKLExpansion = CreateRandomFieldUnique("Short2DKL");
    EXPECT_EQ(twoDimKLExpansion->Name(), "TwoDimKLExpansion");
    
    auto complexTwoDimKLExpansion = CreateRandomFieldUnique("Long2DKL");
    EXPECT_EQ(complexTwoDimKLExpansion->Name(), "ComplexTwoDimKLExpansion");
}

int main(int argc, char **argv) {
  return MppTest(
      MppTestBuilder(argc, argv).
          WithConfigEntry("ConfigVerbose", "0").
          WithConfigEntry("NewtonVerbose", "0").
          WithConfigEntry("LinearVerbose", "0").
          WithConfigEntry("MeshVerbose", "0").
          WithConfigEntry("AssembleVerbose", "0").
          WithRandomInitialized().
          WithPPM()
  ).RUN_ALL_MPP_TESTS();
}