#include "Row.hpp"
#include "TestEnvironment.hpp"

#include <algorithm>

struct RowsParam {
  std::vector<Point> points;
  std::vector<int> dofs;
};

class RowsTest : public Rows, public TestWithParam<RowsParam> {
protected:
  std::vector<Point> points{};
  std::vector<int> dofs{};
  int totalDofs = 0;
  int totalEntries = 0;

  RowsTest() {
    points = GetParam().points;
    std::sort(points.begin(), points.end(),
              [](const Point &p0, const Point &p1) { return p0 < p1; });
    dofs = std::vector<int>(points.size());
    for (int i = 0; i < points.size(); ++i)
      dofs[i] = GetParam().dofs[std::distance(GetParam().points.begin(),
                                              std::find(GetParam().points.begin(),
                                                        GetParam().points.end(), points[i]))];

    mpp_geometry::SetTolerance(1e-14);
    if (dofs.size() != points.size()) dofs = std::vector<int>(points.size(), 1);

    for (int i = 0; i < dofs.size(); ++i)
      totalDofs += dofs[i];
    totalDofs *= totalDofs;
    totalEntries = (points.size() * (points.size() + 1)) / 2;

    for (int i = 0; i < points.size(); ++i)
      Rows::Insert(points[i], dofs[i]);
    for (int i = 1; i < points.size(); ++i)
      for (int j = 0; j < i; ++j)
        AddEntry(points[i], points[j]);
    InitIndices(false);
  }
};

TEST_P(RowsTest, NumberOfRowsTest) { EXPECT_EQ(Rows::size(), points.size()); }

TEST_P(RowsTest, RowNumberOfEntriesTest) {
  for (int i = 0; i < points.size(); ++i) {
    auto r_i = find(points[i]);
    EXPECT_EQ(r_i->second.NumberOfEntries(), i);
    EXPECT_EQ(RowNumberOfEntries(points[i]), i);
  }
}

TEST_P(RowsTest, NumberOfEntriesTotalTest) { EXPECT_EQ(NumberOfEntriesTotal(), totalEntries); }

TEST_P(RowsTest, NumberOfDofsTotalTest) { EXPECT_EQ(NumberOfDofsTotal(), totalDofs); }

TEST_P(RowsTest, RowIdTest) {
  for (int i = 0; i < points.size(); ++i) {
    EXPECT_EQ(RowId(points[i]), i);
    EXPECT_EQ(RowIdChecked(points[i]), i);
  }
}

TEST_P(RowsTest, IndexTest) {
  int k = 0;
  for (int i = 0; i < points.size(); ++i) {
    auto r_i = find(points[i]);
    EXPECT_EQ(Index(r_i, r_i), k * k);

    std::vector<int> idxEntries(i);
    for (int j = 0; j < i; ++j)
      idxEntries[j] = std::distance(r_i->second.EntriesBegin(), r_i->second.FindEntry(points[j]));

    for (int j = 0; j < i; ++j) {
      auto r_j = find(points[j]);
      int m = k * k + dofs[i] * dofs[i];
      for (int s = 0; s < i; ++s)
        if (idxEntries[s] < idxEntries[j]) m += 2 * dofs[i] * dofs[s];

      EXPECT_EQ(Index(r_i, r_j), m);
      EXPECT_EQ(Index(r_j, r_i), m + dofs[i] * dofs[j]);
    }
    k += dofs[i];
  }
}

#if SpaceDimension == 1
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(RowsTest);
#endif

#if SpaceDimension >= 2
INSTANTIATE_TEST_SUITE_P(
    RowsTest, RowsTest,
    Values(RowsParam{{Point(0, 0), Point(1, 0), Point(0, 1), Point(1, 1)}, {1, 2, 3, 4}},
           RowsParam{{Point(0, 1), Point(1, 1), Point(1, 0), Point(0, 2), Point(0, 0)},
                     {5, 1, 2, 3, 2}}));
#endif


int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM();
  return mppTest.RUN_ALL_MPP_TESTS();
}
