#include "Tensor.hpp"
#include "TestEnvironment.hpp"

#include <vector>

template<int sDim>
class IATensorTest : public ::testing::Test {
protected:
  using IATENSOR = TensorT<IAInterval, sDim>;
  using IASYMTENSOR = SymTensorT<IAInterval, sDim>;
  using IAVECTORFIELD = VectorFieldT<IAInterval, sDim>;

  int dim = -1;
  IATENSOR A, B;
  IASYMTENSOR Sym;
  IAVECTORFIELD v;
  IAInterval value;
  double value_d;
  int value_i;

  std::vector<std::vector<IAInterval>> values;
  std::vector<IATENSOR> tensors;
  std::vector<IATENSOR> vfTensors;

  IATensorTest() {
    mpp_geometry::SetTolerance(1e-15);
    dim = A.Dim();
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        A[i][j] = RandomIAInterval();
        B[i][j] = RandomIAInterval();
      }
    for (int i = 0; i < Sym.size(); ++i)
      Sym[i] = RandomIAInterval();
    for (int i = 0; i < dim; ++i)
      v[i] = RandomIAInterval();
    value = RandomNonZeroIAInterval();
    value_d = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();

    values.resize(dim);
    tensors.resize(dim);
    vfTensors.resize(dim);
    for (int i = 0; i < dim; ++i) {
      values[i].resize(dim);
      for (int j = 0; j < dim; ++j)
        values[i][j] = RandomIAInterval();
    }

    vfTensors[0] = IATENSOR(v);
    if constexpr (sDim >= 2) {
      tensors[1] = IATENSOR(values[0][0], values[0][1], values[1][0], values[1][1]);
      vfTensors[1] = IATENSOR(v, v);
    }
    if constexpr (sDim >= 3) {
      tensors[2] = IATENSOR(values[0][0], values[0][1], values[0][2], values[1][0], values[1][1],
                            values[1][2], values[2][0], values[2][1], values[2][2]);
      vfTensors[2] = IATENSOR(v, v, v);
    }
  }

  void checkConstructor() {
    IATENSOR z0;
    IATENSOR z1(A);
    IATENSOR z2(mid(A));
    IATENSOR z3(value);
    IATENSOR z4(value_d);
    IATENSOR z5(value_i);
    IATENSOR z6(Sym);
    IATENSOR z7(mid(Sym));
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(z0[i][j], IAInterval());
        EXPECT_EQ(z1[i][j], A[i][j]);
        EXPECT_EQ(z2[i][j], IAInterval(mid(A[i][j])));
        EXPECT_EQ(z3[i][j], value);
        EXPECT_EQ(z4[i][j], IAInterval(value_d));
        EXPECT_EQ(z5[i][j], IAInterval(value_i));
        EXPECT_EQ(z6[i][j], Sym(i, j));
        EXPECT_EQ(z7[i][j], IAInterval(mid(Sym(i, j))));
      }

    for (int i = 1; i < dim; ++i) {
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          if (k <= i && l <= i) EXPECT_EQ(tensors[i][k][l], values[k][l]);
          else if (k == l) EXPECT_EQ(tensors[i][k][l], IAInterval(1.0));
          else EXPECT_EQ(tensors[i][k][l], IAInterval());
        }
      }
    }

    for (int i = 0; i < dim; ++i) {
      IATENSOR T(i, v);
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          if (k == i) EXPECT_EQ(T[k][l], v[l]);
          else EXPECT_EQ(T[k][l], IAInterval());

          if (k <= i) EXPECT_EQ(vfTensors[i][k][l], v[l]);
          else EXPECT_EQ(vfTensors[i][k][l], IAInterval());
        }
      }
      for (int j = 0; j < dim; ++j) {
        IATENSOR T(i, j, value);
        for (int k = 0; k < dim; ++k) {
          for (int l = 0; l < dim; ++l) {
            if (k == i && l == j) EXPECT_EQ(T[k][l], value);
            else EXPECT_EQ(T[k][l], IAInterval());
          }
        }
      }
    }
  }

  virtual void checkAssignment() {
    IATENSOR z1 = A;
    IATENSOR z2 = mid(A);
    IATENSOR z3(A);
    z3 = value;
    IATENSOR z4(A);
    z4 = value_d;
    IATENSOR z5(A);
    z5 = value_i;
    IATENSOR z6 = Sym;
    IATENSOR z7 = mid(Sym);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(z1[i][j], A[i][j]);
        EXPECT_EQ(z2[i][j], IAInterval(mid(A[i][j])));
        EXPECT_EQ(z3[i][j], value);
        EXPECT_EQ(z4[i][j], IAInterval(value_d));
        EXPECT_EQ(z5[i][j], IAInterval(value_i));
        EXPECT_EQ(z6[i][j], Sym(i, j));
        EXPECT_EQ(z7[i][j], IAInterval(mid(Sym(i, j))));
      }
  }

  virtual void checkSummation() {
    auto r0 = A + B;
    auto r1 = mid(A) + B;
    auto r2 = A + mid(B);
    auto r3 = A + Sym;
    auto r4 = Sym + B;
    auto r5 = A + mid(Sym);
    auto r6 = mid(Sym) + B;
    auto r7 = mid(A) + Sym;
    auto r8 = Sym + mid(B);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r0[i][j], A[i][j] + B[i][j]);
        EXPECT_EQ(r1[i][j], mid(A[i][j]) + B[i][j]);
        EXPECT_EQ(r2[i][j], A[i][j] + mid(B[i][j]));
        EXPECT_EQ(r3[i][j], A[i][j] + Sym(i, j));
        EXPECT_EQ(r4[i][j], Sym(i, j) + B[i][j]);
        EXPECT_EQ(r5[i][j], A[i][j] + mid(Sym(i, j)));
        EXPECT_EQ(r6[i][j], mid(Sym(i, j)) + B[i][j]);
        EXPECT_EQ(r7[i][j], mid(A[i][j]) + Sym(i, j));
        EXPECT_EQ(r8[i][j], Sym(i, j) + mid(B[i][j]));
      }
  }

  virtual void checkDifference() {
    auto r0 = A - B;
    auto r1 = mid(A) - B;
    auto r2 = A - mid(B);
    auto r3 = A - Sym;
    auto r4 = Sym - B;
    auto r5 = A - mid(Sym);
    auto r6 = mid(Sym) - B;
    auto r7 = mid(A) - Sym;
    auto r8 = Sym - mid(B);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r0[i][j], A[i][j] - B[i][j]);
        EXPECT_EQ(r1[i][j], mid(A[i][j]) - B[i][j]);
        EXPECT_EQ(r2[i][j], A[i][j] - mid(B[i][j]));
        EXPECT_EQ(r3[i][j], A[i][j] - Sym(i, j));
        EXPECT_EQ(r4[i][j], Sym(i, j) - B[i][j]);
        EXPECT_EQ(r5[i][j], A[i][j] - mid(Sym(i, j)));
        EXPECT_EQ(r6[i][j], mid(Sym(i, j)) - B[i][j]);
        EXPECT_EQ(r7[i][j], mid(A[i][j]) - Sym(i, j));
        EXPECT_EQ(r8[i][j], Sym(i, j) - mid(B[i][j]));
      }
  }

  virtual void checkAdditiveInverse() {
    auto r = -A;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_EQ(r[i][j], -A[i][j]);
  }

  virtual void checkScalarMultiplication() {
    auto r1 = A * value;
    auto r2 = value * A;
    auto r3 = A * value_d;
    auto r4 = value_d * A;
    auto r5 = A * value_i;
    auto r6 = value_i * A;
    auto r7 = value * mid(A);
    auto r8 = mid(A) * value;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r1[i][j], A[i][j] * value);
        EXPECT_EQ(r2[i][j], A[i][j] * value);
        EXPECT_EQ(r3[i][j], A[i][j] * value_d);
        EXPECT_EQ(r4[i][j], A[i][j] * value_d);
        EXPECT_EQ(r5[i][j], A[i][j] * value_i);
        EXPECT_EQ(r6[i][j], A[i][j] * value_i);
        EXPECT_EQ(r7[i][j], mid(A[i][j]) * value);
        EXPECT_EQ(r8[i][j], mid(A[i][j]) * value);
      }
  }

  virtual void checkScalarDivision() {
    auto r1 = A / value;
    auto r2 = A / value_d;
    auto r3 = A / value_i;
    auto r4 = mid(A) / value;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_EQ(r1[i][j], A[i][j] / value);
        EXPECT_EQ(r2[i][j], A[i][j] / value_d);
        EXPECT_EQ(r3[i][j], A[i][j] / value_i);
        EXPECT_EQ(r4[i][j], mid(A[i][j]) / value);
      }
  }

  virtual void checkMatrixVectorProduct() {
    IAVECTORFIELD r0, r1, r2;
    for (int i = 0; i < dim; ++i)
      for (int k = 0; k < dim; ++k) {
        r0[i] += A[i][k] * v[k];
        r1[i] += A[i][k] * mid(v[k]);
        r2[i] += mid(A[i][k]) * v[k];
      }
    EXPECT_EQ(r0, A * v);
    EXPECT_EQ(r1, A * mid(v));
    EXPECT_EQ(r2, mid(A) * v);
  }

  virtual void checkMatrixMatrixProduct() {
    IATENSOR r0, r1, r2;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k) {
          r0[i][j] += A[i][k] * B[k][j];
          r1[i][j] += A[i][k] * mid(B[k][j]);
          r2[i][j] += mid(A[i][k]) * B[k][j];
        }
    EXPECT_EQ(r0, A * B);
    EXPECT_EQ(r1, A * mid(B));
    EXPECT_EQ(r2, mid(A) * B);
  }

  virtual void checkVectorVectorProduct() {
    IAVECTORFIELD w;
    for (int i = 0; i < w.Dim(); ++i)
      w[i] = RandomIAInterval();
    IATENSOR C0, C1, C2;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        C0[i][j] = v[i] * w[j];
        C1[i][j] = v[i] * mid(w[j]);
        C2[i][j] = mid(v[i]) * w[j];
      }
    EXPECT_EQ(C0, Product(v, w));
    EXPECT_EQ(C1, Product(v, mid(w)));
    EXPECT_EQ(C2, Product(mid(v), w));
  }

  virtual void checkDet() {
    if (A.Dim() == 1) {
      IAInterval d = A[0][0];
      EXPECT_EQ(d, det(A));
    } else if (A.Dim() == 2) {
      IAInterval d = A[0][0] * A[1][1] - A[1][0] * A[0][1];
      EXPECT_EQ(d, det(A));
    } else if (A.Dim() == 3) {
      IAInterval d = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
                     + A[1][0] * (A[2][1] * A[0][2] - A[0][1] * A[2][2])
                     + A[2][0] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);
      EXPECT_EQ(d, det(A));
    } else {
      EXPECT_ANY_THROW(det(A));
    }
  }

  virtual void checkTrace() {
    IAInterval t;
    for (int i = 0; i < dim; ++i)
      t += A[i][i];
    EXPECT_EQ(t, trace(A));
  }

  virtual void checkTranspose() {
    auto r = transpose(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_EQ(r[i][j], A[j][i]);
  }

  virtual void checkFrobeniusProduct() {
    IAInterval sum;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        sum += A[i][j] * B[i][j];
    EXPECT_IAINTERVAL_NEAR(Frobenius(A, B), sum, 1e-14);
  }

  virtual void checkFrobeniusNorm() {
    IAInterval sum;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        sum += sqr(A[i][j]);
    EXPECT_IAINTERVAL_NEAR(norm(A), sqrt(sum), 1e-14);
  }

  virtual void checkDiag() {
    auto r = diag(v);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i == j) EXPECT_EQ(r[i][j], v[i]);
        else EXPECT_EQ(r[i][j], IAInterval(0.0));
    auto d = diag(A);
    for (int i = 0; i < dim; ++i)
      EXPECT_EQ(d[i], A[i][i]);
  }

  virtual void checkSym() {
    IATENSOR C;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i == j) C[i][j] = A[i][j];
        else C[i][j] = 0.5 * (A[i][j] + A[j][i]);
    EXPECT_EQ(C, sym(A));
  }

  virtual void checkSkew() {
    IATENSOR C;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i != j) C[i][j] = 0.5 * (A[i][j] - A[j][i]);
    EXPECT_EQ(C, skew(A));
  }

  virtual void checkEQ() {
    EXPECT_EQ(A, A);
    EXPECT_EQ(B, B);
  }

  virtual void checkNonEQ() {
    EXPECT_NE(A, B);
    EXPECT_NE(B, A);
  }

  virtual void checkSaveLoad() {
    Saver saver("SaveLoadTest");
    saver << A;
    saver.close();
    Loader loader("SaveLoadTest");
    IATENSOR C;
    loader >> C;
    loader.close();
    EXPECT_EQ(A, C);
  }

  virtual void checkMid() {
    auto r = mid(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r[i][j], mid(A[i][j]));
  }

  virtual void checkInf() {
    auto r = inf(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r[i][j], inf(A[i][j]));
  }

  virtual void checkSup() {
    auto r = sup(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r[i][j], sup(A[i][j]));
  }
};

/// To avoid the same code in each test
#define TENSOR_TESTS(tensorClass)                                                                  \
                                                                                                   \
  TEST_F(tensorClass, ConstructorTest) { checkConstructor(); }                                     \
                                                                                                   \
  TEST_F(tensorClass, AssignmentTest) { checkAssignment(); }                                       \
                                                                                                   \
  TEST_F(tensorClass, SummationTest) { checkSummation(); }                                         \
                                                                                                   \
  TEST_F(tensorClass, DifferenceTest) { checkDifference(); }                                       \
                                                                                                   \
  TEST_F(tensorClass, AdditiveInverseTest) { checkAdditiveInverse(); }                             \
                                                                                                   \
  TEST_F(tensorClass, ScalarMultiplicationTest) { checkScalarMultiplication(); }                   \
                                                                                                   \
  TEST_F(tensorClass, ScalarDivisionTest) { checkScalarDivision(); }                               \
                                                                                                   \
  TEST_F(tensorClass, MatrixVectorProductTest) { checkMatrixVectorProduct(); }                     \
                                                                                                   \
  TEST_F(tensorClass, MatrixMatrixProductTest) { checkMatrixMatrixProduct(); }                     \
                                                                                                   \
  TEST_F(tensorClass, VectorVectorProductTest) { checkVectorVectorProduct(); }                     \
                                                                                                   \
  TEST_F(tensorClass, DeterminantTest) { checkDet(); }                                             \
                                                                                                   \
  TEST_F(tensorClass, TraceTest) { checkTrace(); }                                                 \
                                                                                                   \
  TEST_F(tensorClass, TransposeTest) { checkTranspose(); }                                         \
                                                                                                   \
  TEST_F(tensorClass, FrobeniusProductTest) { checkFrobeniusProduct(); }                           \
                                                                                                   \
  TEST_F(tensorClass, FrobeniusNormTest) { checkFrobeniusNorm(); }                                 \
                                                                                                   \
  TEST_F(tensorClass, DiagonalTest) { checkDiag(); }                                               \
                                                                                                   \
  TEST_F(tensorClass, SymmetricTest) { checkSym(); }                                               \
                                                                                                   \
  TEST_F(tensorClass, SkewSymmetricTest) { checkSkew(); }                                          \
                                                                                                   \
  TEST_F(tensorClass, EqualTest) { checkEQ(); }                                                    \
                                                                                                   \
  TEST_F(tensorClass, NonEqualTest) { checkNonEQ(); }                                              \
                                                                                                   \
  TEST_F(tensorClass, SaveLoadTest) { checkSaveLoad(); }                                           \
                                                                                                   \
  TEST_F(tensorClass, MidTest) { checkMid(); }                                                     \
                                                                                                   \
  TEST_F(tensorClass, InfTest) { checkInf(); }                                                     \
                                                                                                   \
  TEST_F(tensorClass, SupTest) { checkSup(); }

using IATensorTest_3 = IATensorTest<3>;

TENSOR_TESTS(IATensorTest_3)

using IATensorTest_2 = IATensorTest<2>;

TENSOR_TESTS(IATensorTest_2)

using IATensorTest_1 = IATensorTest<1>;

TENSOR_TESTS(IATensorTest_1)

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
