#ifndef TESTCOARSEGEOMETRY_HPP
#define TESTCOARSEGEOMETRY_HPP

#include "CoarseGeometry.hpp"
#include "TestEnvironment.hpp"

void EXPECT_GEO_EQ(const CoarseGeometry &geo1, const CoarseGeometry &geo2) {
  auto coord1 = geo1.GetCoordinates();
  auto coord2 = geo2.GetCoordinates();

  EXPECT_EQ(coord1.size(), coord2.size());
  for (const auto &coordinate : coord1) {
    EXPECT_TRUE(std::find(coord2.begin(), coord2.end(), coordinate) != coord2.end());
  }

  auto cells1 = geo1.GetCellIds();
  auto cells2 = geo2.GetCellIds();

  EXPECT_EQ(cells1.size(), cells2.size());
  for (const auto &c : cells1) {
    EXPECT_TRUE(std::find(cells2.begin(), cells2.end(), c) != cells2.end());
  }

  auto faces1 = geo1.GetFaceIds();
  auto faces2 = geo2.GetFaceIds();

  EXPECT_EQ(faces1.size(), faces2.size());
  for (const auto &f : faces1) {
    EXPECT_TRUE(std::find(faces2.begin(), faces2.end(), f) != faces2.end());
  }
}

void EXPECT_DATAGEO_EQ(const CoarseGeometry &geo1, const CoarseGeometry &geo2) {
  EXPECT_GEO_EQ(geo1, geo2);

  auto data1 = geo1.GetVertexDataList();
  auto data2 = geo2.GetVertexDataList();

  EXPECT_EQ(data1.size(), data2.size());
  for (const auto &d : data1) {
    EXPECT_TRUE(std::find(data2.begin(), data2.end(), d) != data2.end());
  }

  data1 = geo1.GetCellDataList();
  data2 = geo2.GetCellDataList();

  EXPECT_EQ(data1.size(), data2.size());
  for (const auto &d : data1) {
    EXPECT_TRUE(std::find(data2.begin(), data2.end(), d) != data2.end());
  }
}

#endif // TESTCOARSEGEOMETRY_HPP
