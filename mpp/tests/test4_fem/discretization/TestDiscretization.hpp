#ifndef TESTDISCRETIZATION_HPP
#define TESTDISCRETIZATION_HPP

#include "IDiscretization.hpp"
#include "MeshesCreator.hpp"
#include "MixedDoF.hpp"
#include "TestEnvironment.hpp"

#include <map>

template<typename T>
class DiscretizationTestT : public TestWithParam<T> {
protected:
  Meshes *meshes = nullptr;
  IDiscretization *disc = nullptr;
  std::string discName = "";
  std::string dofName = "";
  std::map<CELLTYPE, std::map<int, vector<std::string>>> shapeNames = {};
  std::map<CELLTYPE, std::map<int, vector<std::string>>> quadNames = {};
  std::map<CELLTYPE, std::map<int, std::string>> faceQuadNames = {};

  DiscretizationTestT(std::string meshesName, std::string discName, std::string dofName) :
      discName(discName), dofName(dofName) {
    meshes = MeshesCreator(meshesName).Create();
  }

  ~DiscretizationTestT() override {
    if (disc) delete disc;
    if (meshes) delete meshes;
  }

  virtual void checkDiscName() { EXPECT_STREQ(disc->DiscName().c_str(), discName.c_str()); }

  virtual void checkDoFName() {
    EXPECT_STREQ((*disc)({0, 0}).GetDoF().Name().c_str(), dofName.c_str());
  }

  virtual void checkNodalPoints() {
    for (cell c = meshes->fine().cells(); c != meshes->fine().cells_end(); ++c) {
      const auto &shapeNames_c = shapeNames.at(c.ReferenceType());
      for (auto const &[n, names] : shapeNames_c) {
        const IDoF &dof = (*disc)({0, 0}).GetDoF();
        vector<Point> nodalPointsDoF;
        if (typeid(dof) == typeid(MixedDoF)) {
          nodalPointsDoF = MixedDoF::Cast(dof).GetNodalPoints(n, *c);
        } else {
          nodalPointsDoF = dof.GetNodalPoints(*c);
        }
        vector<Point> nodalPointsShape;
        disc->GetShape(*c, n).NodalPoints(*c, nodalPointsShape);
        for (int i = 0; i < nodalPointsDoF.size(); ++i)
          EXPECT_EQ(nodalPointsDoF[i], nodalPointsShape[i]);
      }
    }
  }

  virtual void checkShapeNames() {
    for (cell c = meshes->fine().cells(); c != meshes->fine().cells_end(); ++c) {
      const auto &shapeNames_c = shapeNames.at(c.ReferenceType());
      for (auto const &[n, names] : shapeNames_c) {
        for (int degree = 0; degree < names.size(); ++degree) {
          EXPECT_EQ(disc->GetShape(*c, n, degree).Name(), names[degree]);
          EXPECT_EQ(disc->GetShape(c.ReferenceType(), n, degree).Name(), names[degree]);
        }
      }
    }
  }

  virtual void checkQuadNames() {
    for (cell c = meshes->fine().cells(); c != meshes->fine().cells_end(); ++c) {
      const auto &quadNames_c = quadNames.at(c.ReferenceType());
      for (auto const &[n, names] : quadNames_c) {
        for (int degree = 0; degree < names.size(); ++degree) {
          EXPECT_STREQ(disc->GetQuad(*c, n, degree).Name().c_str(), names[degree].c_str());
          EXPECT_STREQ(disc->GetQuad(c.ReferenceType(), n, degree).Name().c_str(),
                       names[degree].c_str());
        }
      }
    }
  }

  virtual void checkFaceQuadNames() {
    for (cell c = meshes->fine().cells(); c != meshes->fine().cells_end(); ++c) {
      const auto &faceQuadNames_c = faceQuadNames.at(c.ReferenceType());
      for (auto const &[n, name] : faceQuadNames_c) {
        EXPECT_STREQ(disc->GetFaceQuad(*c, n).Name().c_str(), name.c_str());
        EXPECT_STREQ(disc->GetFaceQuad(c.ReferenceType(), n).Name().c_str(), name.c_str());
      }
    }
  }
};

class DiscretizationTest : public DiscretizationTestT<int> {
protected:
  DiscretizationTest(std::string meshesName, std::string discName, std::string dofName) :
      DiscretizationTestT<int>(meshesName, discName, dofName) {}
};

class MixedDiscretizationTest : public DiscretizationTestT<std::pair<int, int>> {
protected:
  MixedDiscretizationTest(std::string meshesName, std::string discName, std::string dofName) :
      DiscretizationTestT<std::pair<int, int>>(meshesName, discName, dofName) {}
};

/// To avoid the same code in each test
#define DISCRETIZATION_TESTS(discTestClass)                                                        \
  TEST_P(discTestClass, DiscNameTest) { checkDiscName(); }                                         \
                                                                                                   \
  TEST_P(discTestClass, DoFNameTest) { checkDoFName(); }                                           \
                                                                                                   \
  TEST_P(discTestClass, NodalPointsTest) { checkNodalPoints(); }                                   \
                                                                                                   \
  TEST_P(discTestClass, ShapeNameTest) { checkShapeNames(); }                                      \
                                                                                                   \
  TEST_P(discTestClass, QuadratureNameTest) { checkQuadNames(); }                                  \
                                                                                                   \
  TEST_P(discTestClass, FaceQuadratureNameTest) { checkFaceQuadNames(); }


#endif // TESTDISCRETIZATION_HPP
