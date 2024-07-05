#include "Tensor.hpp"
#include "TestEnvironment.hpp"

template<int sDim>
class TensorTest : public ::testing::Test {
protected:
  using ROW = TensorRowT<double, sDim>;
  using TENSOR = TensorT<double, sDim>;
  using SYMTENSOR = SymTensorT<double, sDim>;
  using VECTORFIELD = VectorFieldT<double, sDim>;

  int dim = -1;
  TENSOR A, B;
  SYMTENSOR Sym;
  VECTORFIELD v;
  double value;
  int value_i;

  std::vector<std::vector<double>> values;
  std::vector<TENSOR> tensors;
  std::vector<TENSOR> vfTensors;

  TensorTest() {
    mpp_geometry::SetTolerance(1e-15);
    dim = A.Dim();
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        A[i][j] = RandomDouble();
        B[i][j] = RandomDouble();
      }
    for (int i = 0; i < Sym.size(); ++i)
      Sym[i] = RandomDouble();
    for (int i = 0; i < dim; ++i)
      v[i] = RandomDouble();
    value = RandomNonZeroDouble();
    value_i = RandomNonZeroInt();

    values.resize(dim);
    tensors.resize(dim);
    vfTensors.resize(dim);
    for (int i = 0; i < dim; ++i) {
      values[i].resize(dim);
      for (int j = 0; j < dim; ++j)
        values[i][j] = RandomDouble();
    }

    vfTensors[0] = TENSOR(v);
    if constexpr (sDim >= 2) {
      tensors[1] = TENSOR(values[0][0], values[0][1], values[1][0], values[1][1]);
      vfTensors[1] = TENSOR(v, v);
    }
    if constexpr (sDim >= 3) {
      tensors[2] = TENSOR(values[0][0], values[0][1], values[0][2], values[1][0], values[1][1],
                          values[1][2], values[2][0], values[2][1], values[2][2]);
      vfTensors[2] = TENSOR(v, v, v);
    }
  }

  void checkConstructor() {
    TENSOR z0;
    TENSOR z1(A);
    TENSOR z3(value);
    TENSOR z5(value_i);
    TENSOR z6(Sym);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(z0[i][j], 0.0);
        EXPECT_DOUBLE_EQ(z1[i][j], A[i][j]);
        EXPECT_DOUBLE_EQ(z3[i][j], value);
        EXPECT_DOUBLE_EQ(z5[i][j], double(value_i));
        EXPECT_DOUBLE_EQ(z6[i][j], Sym(i, j));
      }

    for (int i = 1; i < dim; ++i) {
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          if (k <= i && l <= i) EXPECT_DOUBLE_EQ(tensors[i][k][l], values[k][l]);
          else if (k == l) EXPECT_DOUBLE_EQ(tensors[i][k][l], 1.0);
          else EXPECT_DOUBLE_EQ(tensors[i][k][l], 0.0);
        }
      }
    }

    for (int i = 0; i < dim; ++i) {
      TENSOR T(i, v);
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          if (k == i) EXPECT_DOUBLE_EQ(T[k][l], v[l]);
          else EXPECT_DOUBLE_EQ(T[k][l], 0.0);

          if (k <= i) EXPECT_DOUBLE_EQ(vfTensors[i][k][l], v[l]);
          else EXPECT_DOUBLE_EQ(vfTensors[i][k][l], 0.0);
        }
      }
      for (int j = 0; j < dim; ++j) {
        TENSOR T(i, j, value);
        for (int k = 0; k < dim; ++k) {
          for (int l = 0; l < dim; ++l) {
            if (k == i && l == j) EXPECT_DOUBLE_EQ(T[k][l], value);
            else EXPECT_DOUBLE_EQ(T[k][l], 0.0);
          }
        }
      }
    }
  }

  virtual void checkAssignment() {
    TENSOR z1 = A;
    TENSOR z3(A);
    z3 = value;
    TENSOR z5(A);
    z5 = value_i;
    TENSOR z6 = Sym;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(z1[i][j], A[i][j]);
        EXPECT_DOUBLE_EQ(z3[i][j], value);
        EXPECT_DOUBLE_EQ(z5[i][j], double(value_i));
        EXPECT_DOUBLE_EQ(z6[i][j], Sym(i, j));
      }
  }

  virtual void checkSummation() {
    auto r0 = A + B;
    auto r3 = A + Sym;
    auto r4 = Sym + B;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r0[i][j], A[i][j] + B[i][j]);
        EXPECT_DOUBLE_EQ(r3[i][j], A[i][j] + Sym(i, j));
        EXPECT_DOUBLE_EQ(r4[i][j], Sym(i, j) + B[i][j]);
      }
  }

  virtual void checkDifference() {
    auto r0 = A - B;
    auto r3 = A - Sym;
    auto r4 = Sym - B;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r0[i][j], A[i][j] - B[i][j]);
        EXPECT_DOUBLE_EQ(r3[i][j], A[i][j] - Sym(i, j));
        EXPECT_DOUBLE_EQ(r4[i][j], Sym(i, j) - B[i][j]);
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
    auto r5 = A * value_i;
    auto r6 = value_i * A;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r1[i][j], A[i][j] * value);
        EXPECT_DOUBLE_EQ(r2[i][j], A[i][j] * value);
        EXPECT_DOUBLE_EQ(r5[i][j], A[i][j] * value_i);
        EXPECT_DOUBLE_EQ(r6[i][j], A[i][j] * value_i);
      }
  }

  virtual void checkScalarDivision() {
    auto r1 = A / value;
    auto r3 = A / value_i;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        EXPECT_DOUBLE_EQ(r1[i][j], A[i][j] / value);
        EXPECT_DOUBLE_EQ(r3[i][j], A[i][j] / value_i);
      }
  }

  virtual void checkMatrixVectorProduct() {
    VECTORFIELD r0;
    for (int i = 0; i < dim; ++i)
      for (int k = 0; k < dim; ++k) {
        r0[i] += A[i][k] * v[k];
      }
    EXPECT_EQ(r0, A * v);
  }

  virtual void checkMatrixMatrixProduct() {
    TENSOR r0;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k) {
          r0[i][j] += A[i][k] * B[k][j];
        }
    EXPECT_EQ(r0, A * B);
  }

  virtual void checkRowMatrixProduct() {
    for (int k = 0; k < dim; ++k) {
      ROW r(k, v);
      TENSOR R(k, v);
      EXPECT_EQ(TENSOR(r * A), R * A);
    }
  }

  virtual void checkMatrixRowProduct() {
    for (int k = 0; k < dim; ++k) {
      ROW r(k, v);
      TENSOR R(k, v);
      EXPECT_EQ(A * r, A * R);
    }
  }

  virtual void checkVectorVectorProduct() {
    VECTORFIELD w;
    for (int i = 0; i < w.Dim(); ++i)
      w[i] = RandomDouble();
    TENSOR C0;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        C0[i][j] = v[i] * w[j];
      }
    EXPECT_EQ(C0, Product(v, w));
  }

  virtual void checkDet() {
    if (A.Dim() == 1) {
      double d = A[0][0];
      EXPECT_DOUBLE_EQ(d, det(A));
    } else if (A.Dim() == 2) {
      double d = A[0][0] * A[1][1] - A[1][0] * A[0][1];
      EXPECT_DOUBLE_EQ(d, det(A));
    } else if (A.Dim() == 3) {
      double d = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
                 + A[1][0] * (A[2][1] * A[0][2] - A[0][1] * A[2][2])
                 + A[2][0] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);
      EXPECT_DOUBLE_EQ(d, det(A));
    } else {
      EXPECT_ANY_THROW(det(A));
    }
  }

  virtual void checkTrace() {
    double t = 0.0;
    for (int i = 0; i < dim; ++i)
      t += A[i][i];
    EXPECT_DOUBLE_EQ(t, trace(A));
  }

  virtual void checkTranspose() {
    auto r = transpose(A);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        EXPECT_DOUBLE_EQ(r[i][j], A[j][i]);
  }

  virtual void checkFrobeniusProduct() {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        sum += A[i][j] * B[i][j];
    EXPECT_DOUBLE_EQ(Frobenius(A, B), sum);
  }

  virtual void checkFrobeniusNorm() {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        sum += A[i][j] * A[i][j];
    EXPECT_DOUBLE_EQ(norm(A), sqrt(sum));
  }

  virtual void checkDiag() {
    auto r = diag(v);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i == j) EXPECT_DOUBLE_EQ(r[i][j], v[i]);
        else EXPECT_DOUBLE_EQ(r[i][j], 0.0);
    auto d = diag(A);
    for (int i = 0; i < dim; ++i)
      EXPECT_DOUBLE_EQ(d[i], A[i][i]);
  }

  virtual void checkSym() {
    TENSOR C;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i == j) C[i][j] = A[i][j];
        else C[i][j] = 0.5 * (A[i][j] + A[j][i]);
    EXPECT_EQ(C, sym(A));
  }

  virtual void checkSkew() {
    TENSOR C;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i != j) C[i][j] = 0.5 * (A[i][j] - A[j][i]);
    EXPECT_EQ(C, skew(A));
  }

  virtual void checkCross() {
    if (A.Dim() == 3) {
      TENSOR C = 0.5 * Cross(A, A);
      TENSOR Cof = det(A) * transposeInvert(A);

      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          EXPECT_NEAR(C[i][j], Cof[i][j], 1e-6);
        };
      }
    }
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
    TENSOR C;
    loader >> C;
    loader.close();
    EXPECT_EQ(A, C);
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
  TEST_F(tensorClass, MatrixRowProductTest) { checkMatrixRowProduct(); }                           \
                                                                                                   \
  TEST_F(tensorClass, RowMatrixProductTest) { checkRowMatrixProduct(); }                           \
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
  TEST_F(tensorClass, CrossProductTest) { checkCross(); }                                          \
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
  TEST_F(tensorClass, SaveLoadTest) { checkSaveLoad(); }

using TensorTest_3 = TensorTest<SpaceDimension>;

TENSOR_TESTS(TensorTest_3)

/*using TensorTest_2 = TensorTest<2>;

TENSOR_TESTS(TensorTest_2)

using TensorTest_1 = TensorTest<1>;

TENSOR_TESTS(TensorTest_1)*/

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
