#ifndef TESTPARALLELIZATION_HPP
#define TESTPARALLELIZATION_HPP

#include "TestEnvironment.hpp"

struct BroadcastStruct {
  int integer;

  BroadcastStruct(int _integer) : integer(_integer) {}

  bool operator==(const BroadcastStruct other) { return (other.integer == integer); }
};

class TestPPM : public Test {
protected:
  int commSplit = 0;

  void SetUp() override { PPM->ClearCommunicators(); }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(0);
  }

  virtual void checkIsInitialized() { EXPECT_TRUE(PPM->IsInitialized()); }

  virtual void checkBroadcastDouble() {
    double a = 1.0;
    double b;
    if (PPM->Master(commSplit)) PPM->Broadcast(&a, sizeof(double), commSplit);
    if (!PPM->Master(commSplit)) {
      PPM->Broadcast(&b, sizeof(double), commSplit);
      EXPECT_DOUBLE_EQ(a, b);
    }
  }

  virtual void checkBroadcastDouble2() {
    double a;
    if (PPM->Master(commSplit)) PPM->BCastOnCommSplit(1.0, commSplit);
    if (!PPM->Master(commSplit)) {
      PPM->BCastOnCommSplit(a, commSplit);
      EXPECT_DOUBLE_EQ(a, 1.0);
    }
  }

  virtual void checkBroadcastDouble3() {
    double a;
    if (PPM->Master(0)) PPM->BroadcastDouble(1.0);
    if (!PPM->Master(0)) {
      a = PPM->BroadcastDouble();
      EXPECT_DOUBLE_EQ(a, 1.0);
    }
  }

  virtual void checkBroadcastInt() {
    int a = 1;
    int b;
    if (PPM->Master(commSplit)) PPM->Broadcast(&a, sizeof(int), commSplit);
    if (!PPM->Master(commSplit)) {
      PPM->Broadcast(&b, sizeof(double), commSplit);
      EXPECT_EQ(a, b);
    }
  }

  virtual void checkBroadcastInt2() {
    int a;
    if (PPM->Master(commSplit)) PPM->BCastOnCommSplit(1, commSplit);
    if (!PPM->Master(commSplit)) {
      PPM->BCastOnCommSplit(a, commSplit);
      EXPECT_EQ(a, 1);
    }
  }

  virtual void checkBroadcastInt3() {
    int a;
    if (PPM->Master(0)) PPM->BroadcastInt(1);
    if (!PPM->Master(0)) {
      a = PPM->BroadcastInt();
      EXPECT_EQ(a, 1);
    }
  }

  virtual void checkBroadcast() {
    BroadcastStruct structToBroadcast(0);
    BroadcastStruct structFromBroadcast(1);
    if (PPM->Master(commSplit))
      PPM->Broadcast(&structToBroadcast, sizeof(BroadcastStruct), commSplit);
    if (!PPM->Master(commSplit)) {
      PPM->Broadcast(&structFromBroadcast, sizeof(BroadcastStruct), commSplit);
      EXPECT_TRUE(structFromBroadcast == structToBroadcast);
    }
  }

  virtual void checkBroadcast2() {
    const BroadcastStruct structToBroadcast(0);
    BroadcastStruct structFromBroadcast(1);
    if (PPM->Master(commSplit)) PPM->BCastOnCommSplit(structToBroadcast, commSplit);
    if (!PPM->Master(commSplit)) {
      PPM->BCastOnCommSplit(structFromBroadcast, commSplit);
      EXPECT_TRUE(structFromBroadcast == structToBroadcast);
    }
  }

  virtual void checkProc() {
    if (PPM->Master(commSplit)) EXPECT_EQ(PPM->Proc(commSplit), 0);
    if (!PPM->Master(commSplit)) EXPECT_NE(PPM->Proc(commSplit), 0);
  }

  virtual void checkSize() {
    int testSize;
    MPI_Comm_size(MPI_COMM_WORLD, &testSize);
    EXPECT_EQ(PPM->Size(commSplit), testSize);
  }

  virtual void checkColor() { EXPECT_EQ(PPM->Color(commSplit), 0); }

  virtual void checkMaxCommSplit() { EXPECT_EQ(commSplit, PPM->MaxCommSplit()); }

  virtual void checkClearCommunicators() {
    PPM->ClearCommunicators();
    EXPECT_EQ(PPM->CommsSize(), 1);
  }

  virtual void checkIntegerSum() {
    int a = 1;
    PPM->Sum(&a, 1, commSplit);
    EXPECT_EQ(a, PPM->Size(commSplit));
  }

  virtual void checkIntegerSum2() {
    int a = 1;
    EXPECT_EQ(PPM->SumOnCommSplit(a, commSplit), PPM->Size(commSplit));
  }

  virtual void checkDoubleSum() {
    double a = 1.0;
    PPM->Sum(&a, 1, commSplit);
    EXPECT_DOUBLE_EQ(a, (double)PPM->Size(commSplit));
  }

  virtual void checkDoubleSum2() {
    double a = 1.0;
    EXPECT_DOUBLE_EQ(PPM->SumOnCommSplit(a, commSplit), (double)PPM->Size(commSplit));
  }

  virtual void checkSizeTypeSum() {
    size_t a = sizeof(int);
    PPM->Sum(&a, 1, commSplit);
    EXPECT_EQ(a, (size_t)sizeof(int) * PPM->Size(commSplit));
  }

  virtual void checkSizeTypeSum2() {
    size_t a = sizeof(int);
    EXPECT_EQ(PPM->SumOnCommSplit(a, commSplit), (size_t)a * PPM->Size(commSplit));
  }

  virtual void checkComplexSum() {
    std::complex<double> a{1.0, 1.0};
    std::complex<double> b{std::real(a) * PPM->Size(commSplit),
                           std::imag(a) * PPM->Size(commSplit)};
    PPM->Sum(&a, 1, commSplit);
    EXPECT_EQ(a, b);
  }

  virtual void checkIntegerMin() {
    int a;
    if (PPM->Master(commSplit)) a = 1;
    if (!PPM->Master(commSplit)) a = 0;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Min(a, commSplit), 0);
    else EXPECT_EQ(PPM->Min(a, commSplit), 1);

    if (PPM->Master(commSplit)) a = 0;
    if (!PPM->Master(commSplit)) a = 1;

    EXPECT_EQ(PPM->Min(a, commSplit), 0);
  }

  virtual void checkDoubleMin() {
    double a;
    if (PPM->Master(commSplit)) a = 1.0;
    if (!PPM->Master(commSplit)) a = 0.0;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Min(a, commSplit), 0.0);
    else EXPECT_EQ(PPM->Min(a, commSplit), 1.0);

    if (PPM->Master(commSplit)) a = 0.0;
    if (!PPM->Master(commSplit)) a = 1.0;

    EXPECT_EQ(PPM->Min(a, commSplit), 0.0);
  }

  virtual void checkUnsignedIntegerMin() {
    long unsigned int a;
    if (PPM->Master(commSplit)) a = 1;
    if (!PPM->Master(commSplit)) a = 0;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Min(a, commSplit), 0);
    else EXPECT_EQ(PPM->Min(a, commSplit), 1);

    if (PPM->Master(commSplit)) a = 0;
    if (!PPM->Master(commSplit)) a = 1;

    EXPECT_EQ(PPM->Min(a, commSplit), 0);
  }

  virtual void checkIntegerMax() {
    int a;
    if (PPM->Master(commSplit)) a = 1;
    if (!PPM->Master(commSplit)) a = 0;

    EXPECT_EQ(PPM->Max(a, commSplit), 1);

    if (PPM->Master(commSplit)) a = 0;
    if (!PPM->Master(commSplit)) a = 1;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Max(a, commSplit), 1);
    else EXPECT_EQ(PPM->Max(a, commSplit), 0);
  }

  virtual void checkDoubleMax() {
    double a;
    if (PPM->Master(commSplit)) a = 1.0;
    if (!PPM->Master(commSplit)) a = 0.0;

    EXPECT_DOUBLE_EQ(PPM->Max(a, commSplit), 1.0);

    if (PPM->Master(commSplit)) a = 0.0;
    if (!PPM->Master(commSplit)) a = 1.0;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Max(a, commSplit), 1.0);
    else EXPECT_EQ(PPM->Max(a, commSplit), 0.0);
  }

  virtual void checkUnsignedIntMax() {
    long unsigned int a;
    if (PPM->Master(commSplit)) a = 1;
    if (!PPM->Master(commSplit)) a = 0;

    EXPECT_EQ(PPM->Max(a, commSplit), 1);

    if (PPM->Master(commSplit)) a = 0;
    if (!PPM->Master(commSplit)) a = 1;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Max(a, commSplit), 1);
    else EXPECT_EQ(PPM->Max(a, commSplit), 0);
  }

  virtual void checkAnd() {
    bool b;
    if (PPM->Master(commSplit)) b = true;
    if (!PPM->Master(commSplit)) b = false;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->And(b, 0), false);
    else EXPECT_EQ(PPM->And(b, 0), true);

    if (PPM->Master(commSplit)) b = false;
    if (!PPM->Master(commSplit)) b = true;

    EXPECT_EQ(PPM->And(b, 0), false);
  }

  virtual void checkOr() {
    bool b;
    if (PPM->Master(commSplit)) b = true;
    if (!PPM->Master(commSplit)) b = false;

    EXPECT_EQ(PPM->Or(b), true);

    if (PPM->Master(commSplit)) b = false;
    if (!PPM->Master(commSplit)) b = true;

    if (PPM->Size(commSplit) != 1) EXPECT_EQ(PPM->Or(b), true);
    else EXPECT_EQ(PPM->Or(b), false);
  }
};

class TestPPMWithSplit : public TestPPM {
protected:
  void SetUp() override {
    PPM->ClearCommunicators();
    commSplit = 1;
    PPM->SplitCommunicators();
  }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(0);
  }

  virtual void checkSize() override {
    if (PPM->Size(0) != 1) EXPECT_EQ(PPM->Size(1), PPM->Size(0) / 2);
    else EXPECT_EQ(PPM->Size(1), PPM->Size(0));
  }

  virtual void checkColor() override {
    EXPECT_EQ(PPM->Proc(1) + PPM->Color(1) * PPM->Size(1), PPM->Proc(0));
  }
};

class TestPPMWithDoubleSplit : public TestPPM {
protected:
  void SetUp() override {
    PPM->ClearCommunicators();
    commSplit = 2;
    PPM->SplitCommunicators();
    PPM->SplitCommunicators();
  }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(0);
  }

  void checkSize() override {
    EXPECT_EQ(PPM->Size(1), PPM->Size(0) / 2);
    EXPECT_EQ(PPM->Size(2), PPM->Size(1) / 2);
    EXPECT_EQ(PPM->Size(2), PPM->Size(0) / 4);
  }

  void checkColor() override {
    EXPECT_EQ(PPM->Proc(2) + PPM->Color(2) * PPM->Size(2), PPM->Proc(0));
  }
};

class TestPPMWithFullSplit : public TestPPM {
protected:
  void SetUp() override {
    PPM->ClearCommunicators();
    PPM->FullSplit();
    commSplit = PPM->MaxCommSplit();
    PPM->Barrier(0);
  }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(0);
  }

  void checkSize() override { EXPECT_EQ(PPM->Size(commSplit), 1); }

  void checkColor() override {
    //        EXPECT_EQ(PPM->Proc(2) + PPM->Color(2) * PPM->Size(2), PPM->Proc(1));
    EXPECT_EQ(PPM->Color(commSplit), PPM->Proc(0));
    //        EXPECT_EQ(PPM->Color(commSplit), 0);
  }
};

/// To avoid the same code in each test
#define PPM_TESTS(PPMTestClass)                                                                    \
                                                                                                   \
  TEST_F(PPMTestClass, TestIsInitialized) { checkIsInitialized(); }                                \
  TEST_F(PPMTestClass, TestBroadcastDouble) {                                                      \
    checkBroadcastDouble();                                                                        \
    checkBroadcastDouble2();                                                                       \
    checkBroadcastDouble3();                                                                       \
  }                                                                                                \
  TEST_F(PPMTestClass, TestBroadcastInt) {                                                         \
    checkBroadcastInt();                                                                           \
    checkBroadcastInt2();                                                                          \
    checkBroadcastInt3();                                                                          \
  }                                                                                                \
  TEST_F(PPMTestClass, TestBroadcast) {                                                            \
    checkBroadcast();                                                                              \
    checkBroadcast2();                                                                             \
  }                                                                                                \
  TEST_F(PPMTestClass, TestProc) { checkProc(); }                                                  \
  TEST_F(PPMTestClass, TestSize) { checkSize(); }                                                  \
  TEST_F(PPMTestClass, TestColor) { checkColor(); }                                                \
  TEST_F(PPMTestClass, TestMaxCommSplit) { checkMaxCommSplit(); }                                  \
  TEST_F(PPMTestClass, TestClearCommunicators) { checkClearCommunicators(); }                      \
  TEST_F(PPMTestClass, TestIntegerSum) {                                                           \
    checkIntegerSum();                                                                             \
    checkIntegerSum2();                                                                            \
  }                                                                                                \
  TEST_F(PPMTestClass, TestDoubleSum) {                                                            \
    checkDoubleSum();                                                                              \
    checkDoubleSum2();                                                                             \
  }                                                                                                \
  TEST_F(PPMTestClass, TestSizeTypeSum) {                                                          \
    checkSizeTypeSum();                                                                            \
    checkSizeTypeSum2();                                                                           \
  }                                                                                                \
  TEST_F(PPMTestClass, TestComplexSum) { checkComplexSum(); }                                      \
  TEST_F(PPMTestClass, TestIntegerMin) { checkIntegerMin(); }                                      \
  TEST_F(PPMTestClass, TestDoubleMin) { checkDoubleMin(); }                                        \
  TEST_F(PPMTestClass, TestUnsignedIntegerMin) { checkUnsignedIntegerMin(); }                      \
  TEST_F(PPMTestClass, TestIntegerMax) { checkIntegerMax(); }                                      \
  TEST_F(PPMTestClass, TestDoubleMax) { checkDoubleMax(); }                                        \
  TEST_F(PPMTestClass, TestUnsignedIntMax) { checkUnsignedIntMax(); }                              \
  TEST_F(PPMTestClass, TestAnd) { checkAnd(); }                                                    \
  TEST_F(PPMTestClass, TestOr) { checkOr(); }


#endif // TESTPARALLELIZATION_HPP
