#include "TestCheckOrientation.hpp"
#include "CheckOrientation.hpp"

#include "TestEnvironment.hpp"

using namespace ::testing;

static std::vector<std::vector<int>> admissibleCellsInt{{0, 1}};

static std::vector<std::vector<int>> admissibleCellsTri{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};

static std::vector<std::vector<int>> admissibleCellsQuad{{0, 1, 2, 3},
                                                         {1, 2, 3, 0},
                                                         {2, 3, 0, 1},
                                                         {3, 0, 1, 2}};

static std::vector<std::vector<int>> admissibleCellsTet{{0, 1, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3},
                                                        {1, 3, 2, 0}, {3, 2, 1, 0}, {2, 1, 3, 0},
                                                        {0, 2, 3, 1}, {2, 3, 0, 1}, {3, 0, 2, 1},
                                                        {0, 3, 1, 2}, {3, 1, 0, 2}, {1, 0, 3, 2}};

static std::vector<std::vector<int>>
    admissibleCellsHex{{0, 1, 2, 3, 4, 5, 6, 7}, {0, 3, 7, 4, 1, 2, 6, 5},
                       {0, 4, 5, 1, 3, 7, 6, 2}, {1, 0, 4, 5, 2, 3, 7, 6},
                       {1, 2, 3, 0, 5, 6, 7, 4}, {1, 5, 6, 2, 0, 4, 7, 3},
                       {2, 1, 5, 6, 3, 0, 4, 7}, {2, 3, 0, 1, 6, 7, 4, 5},
                       {2, 6, 7, 3, 1, 5, 4, 0}, {3, 0, 1, 2, 7, 4, 5, 6},
                       {3, 2, 6, 7, 0, 1, 5, 4}, {3, 7, 4, 0, 2, 6, 5, 1},
                       {4, 0, 3, 7, 5, 1, 2, 6}, {4, 5, 1, 0, 7, 6, 2, 3},
                       {4, 7, 6, 5, 0, 3, 2, 1}, {5, 1, 0, 4, 6, 2, 3, 7},
                       {5, 4, 7, 6, 1, 0, 3, 2}, {5, 6, 2, 1, 4, 7, 3, 0},
                       {6, 2, 1, 5, 7, 3, 0, 4}, {6, 5, 4, 7, 2, 1, 0, 3},
                       {6, 7, 3, 2, 5, 4, 0, 1}, {7, 3, 2, 6, 4, 0, 1, 5},
                       {7, 4, 0, 3, 6, 5, 1, 2}, {7, 6, 5, 4, 3, 2, 1, 0}};

const std::vector<std::vector<int>> &getAdmissibleCells(CELLTYPE type) {
  switch (type) {
  case INTERVAL:
    return admissibleCellsInt;
  case TRIANGLE:
    return admissibleCellsTri;
  case QUADRILATERAL:
    return admissibleCellsQuad;
  case TETRAHEDRON:
    return admissibleCellsTet;
  case HEXAHEDRON:
    return admissibleCellsHex;
  default:
    THROW("cell type not implemented")
  }
}

class TestCheckOrientation : public TestWithParam<CheckOrientationTestParameters> {
protected:
  CELLTYPE type;
  std::vector<std::vector<Point>> admissibleCornersList{};
  std::vector<std::vector<Point>> cornersToCheck{};

  TestCheckOrientation(CELLTYPE cellType) : type(cellType) {
    const std::vector<Point> &corners = GetParam().corners;
    const std::vector<std::vector<int>> &admissibleCells = getAdmissibleCells(type);
    for (const auto &admissibleCell : admissibleCells) {
      std::vector<Point> tmp{};
      for (auto id : admissibleCell) {
        tmp.push_back(corners[id]);
      }
      admissibleCornersList.push_back(tmp);
    }
  }

  ~TestCheckOrientation() override {}

  void check() {
    std::vector<Point> corners = GetParam().corners;
    std::cout << "before sort " << corners[0] << corners[1] << endl;
    std::sort(corners.begin(), corners.end());
    std::cout << "after sort " << corners[0] << corners[1] << endl;
    do {
      std::cout << "permut: " << corners[0] << corners[1] << endl;
      std::vector<Point> testCorners = corners;
      std::vector<int> testIndices(corners.size());
      std::iota(testIndices.begin(), testIndices.end(), 0);

      if (GetParam().expectSuccess) {
        CheckOrientation::Check(type, testCorners, testIndices);

        EXPECT_TRUE(checkPoints(testCorners));
        checkIndices(corners, testCorners, testIndices);
      } else {
        EXPECT_THROW(CheckOrientation::Check(type, testCorners, testIndices), MppException);
      }
    } while (std::next_permutation(corners.begin(), corners.end()));
  }

  bool checkPoints(const std::vector<Point> &checkedCorners) {
    for (const auto &admissibleCorners : admissibleCornersList) {
      if (admissibleCorners.size() != checkedCorners.size()) continue;

      bool tmp = true;
      for (int i = 0; i < admissibleCorners.size(); ++i) {
        tmp = tmp && (admissibleCorners[i] == checkedCorners[i]);
      }
      if (tmp) return true;
    }
    return false;
  }

  void checkIndices(const std::vector<Point> &corners, const std::vector<Point> &checkedCorners,
                    const std::vector<int> &testIndices) {
    for (int i = 0; i < corners.size(); ++i) {
      EXPECT_EQ(corners[testIndices[i]], checkedCorners[i]);
    }
  }
};

/// To avoid the same code in each test
#define CHECKORIENTATION_TESTS(testclass, testCases)                                               \
                                                                                                   \
  TEST_P(testclass, TestCheckOrientation) { check(); }                                             \
  INSTANTIATE_TEST_SUITE_P(TestCheckOrientation, testclass, ValuesIn(testCases));

class TestCheckOrientationInterval : public TestCheckOrientation {
public:
  TestCheckOrientationInterval() : TestCheckOrientation(INTERVAL) {}
};

CHECKORIENTATION_TESTS(TestCheckOrientationInterval, TestCasesInt)

#if SpaceDimension >= 2

class TestCheckOrientationTriangle : public TestCheckOrientation {
public:
  TestCheckOrientationTriangle() : TestCheckOrientation(TRIANGLE) {}
};

CHECKORIENTATION_TESTS(TestCheckOrientationTriangle, TestCasesTri)

class TestCheckOrientationQuadrilateral : public TestCheckOrientation {
public:
  TestCheckOrientationQuadrilateral() : TestCheckOrientation(QUADRILATERAL) {}
};

CHECKORIENTATION_TESTS(TestCheckOrientationQuadrilateral, TestCasesQuad)

#endif
#if SpaceDimension >= 3

class TestCheckOrientationTetrahedron : public TestCheckOrientation {
public:
  TestCheckOrientationTetrahedron() : TestCheckOrientation(TETRAHEDRON) {}
};

CHECKORIENTATION_TESTS(TestCheckOrientationTetrahedron, TestCasesTet)

class TestCheckOrientationHexahedron : public TestCheckOrientation {
public:
  TestCheckOrientationHexahedron() : TestCheckOrientation(HEXAHEDRON) {}
};

CHECKORIENTATION_TESTS(TestCheckOrientationHexahedron, TestCasesHex)

#endif

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}