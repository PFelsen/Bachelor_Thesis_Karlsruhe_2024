#include "FFT.hpp"
#include "TestEnvironment.hpp"
#include "m++.hpp"

class FFTTest : public Test {
protected:
  RVector SignalInput = RVector(0, 1000);
  CVector SignalTransformed = CVector(0, 1000);

  FFTTest() {}

  void AddSignalFreq(RVector &SignalInput, int freq, int shift) {
    double dt = 1.0 / SignalInput.size() * 2.0 * Pi * freq;
    for (int i = 0; i < SignalInput.size(); ++i) {
      SignalInput[i] += sin(dt * (i - shift));
    }
  }
};

TEST_F(FFTTest, FFTInversionTest) {
  AddSignalFreq(SignalInput, 10, 0);
  AddSignalFreq(SignalInput, 23, 201);
  AddSignalFreq(SignalInput, 1, 27);
  RVector Tmp(SignalInput);
  FFT::RealToComplexVector(SignalInput, SignalTransformed);
  FFT::InvComplexToRealVector(SignalTransformed, SignalInput);
  EXPECT_NEAR((SignalInput - Tmp).norm(), 0, 1e-14);
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithScreenLogging().WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}