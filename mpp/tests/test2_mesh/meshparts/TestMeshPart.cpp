#include "TestMeshPart.hpp"

TEST(MeshPartTest, InsertionOfNonExisting) {
  auto mp = MeshPart<int, double>();
  mp.Insert(1, 1.0);

  EXPECT_EQ(mp.size(), 1);
  EXPECT_DOUBLE_EQ(mp.Find(1)->second, 1.0);
}

TEST(MeshPartTest, InsertionOfExisting) {
  auto mp = MeshPart<int, double>();
  mp.Insert(1, 1.0);
  mp.Insert(1, 2.0);

  EXPECT_EQ(mp.size(), 1);
  EXPECT_DOUBLE_EQ(mp.Find(1)->second, 1.0);
}

TEST(MeshPartTest, ReplaceNonExisting) {
  auto mp = MeshPart<int, double>();
  mp.Replace(1, 1.0);

  EXPECT_EQ(mp.size(), 1);
  EXPECT_DOUBLE_EQ(mp.Find(1)->second, 1.0);
}

TEST(MeshPartTest, ReplaceExisting) {
  auto mp = MeshPart<int, double>();
  mp.Insert(1, 1.0);
  mp.Replace(1, 2.0);

  EXPECT_EQ(mp.size(), 1);
  EXPECT_DOUBLE_EQ(mp.Find(1)->second, 2.0);
}

TEST(MeshPartTest, DeletionOfExisting) {
  auto mp = MeshPart<int, double>();
  mp.Insert(1, 1.0);
  mp.Remove(1);

  EXPECT_EQ(mp.size(), 0);
}

TEST(MeshPartTest, DeletionOfNonExisting) {
  auto mp = MeshPart<int, double>();
  mp.Insert(1, 1.0);
  mp.Remove(0);

  EXPECT_EQ(mp.size(), 1);
}

TEST(MeshPartTest, ForEachLoop) {
  auto mp = MeshPart<int, double>();
  for (int i = 0; i < 10; ++i) {
    mp.Insert(i, (double)i);
  }
  EXPECT_EQ(mp.size(), 10);

  for (const auto &d : mp) {
    EXPECT_DOUBLE_EQ(double(d.first), d.second);
  }
}

TEST(MeshPartTest, BeginEndLoop) {
  auto mp = MeshPart<int, double>();
  for (int i = 0; i < 10; ++i) {
    mp.Insert(i, (double)i);
  }
  EXPECT_EQ(mp.size(), 10);

  for (auto d = mp.Begin(); d != mp.End(); ++d) {
    EXPECT_DOUBLE_EQ(double(d->first), d->second);
  }
}

int main(int argc, char **argv) {
  // InitGoogleTest(&argc, argv);
  // RUN_ALL_TESTS();
  auto mpptest = MppTest(MppTestBuilder(argc, argv));
  mpptest.RUN_ALL_MPP_TESTS();
}