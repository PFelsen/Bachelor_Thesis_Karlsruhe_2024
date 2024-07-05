#ifndef SPACETIME_TESTREFINEMESH_HPP
#define SPACETIME_TESTREFINEMESH_HPP

#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"

struct TestConfig {
  std::string meshName;
  std::string meshNameRefinied;
  std::string distName;
  std::string olapName;
};

class TestRefineMesh : public TestWithParam<TestConfig> {
protected:
  std::string meshName;
  std::string meshNameRefinied;
  std::string distName;
  std::unique_ptr<Meshes> meshes0;
  std::unique_ptr<Meshes> meshes1;

  TestRefineMesh() :
      meshName(GetParam().meshName), meshNameRefinied(GetParam().meshNameRefinied),
      distName(GetParam().distName) {

    meshes0 =
        MeshesCreator(meshName).WithDistribute(distName).WithPLevel(0).WithLevel(0).CreateUnique();

    meshes1 = MeshesCreator(meshNameRefinied)
                  .WithDistribute(distName)
                  .WithPLevel(0)
                  .WithLevel(0)
                  .CreateUnique();
  }

  bool areEqual(const Mesh &m1, const Mesh &m2) {
    std::map<string, bool> equalities;
    equalities["procSets"] = m1.GetProcSets() == m2.GetProcSets();
    equalities["identifySets"] = m1.identifySets == m2.identifySets;

    equalities["cells"] = true;
    for (cell c = m1.cells(); c != m1.cells_end(); c++) {
      cell c2 = m2.find_cell(c());
      if (c2 == m2.cells_end()) {
        equalities["cells"] = false;
      } else {
        mout << c.copyCorners() << endl;
        mout << c2.copyCorners() << endl;
        equalities["cells"] &= c.copyCorners() == c2.copyCorners();
      }
    }
    equalities["cellcount"] = m1.CellCount() == m2.CellCount();
    equalities["edges"] = true;
    for (edge e = m1.edges(); e != m1.edges_end(); e++) {
      edge e2 = m2.find_edge(e());
      if (e2 == m2.edges_end()) {
        equalities["edges"] = false;
      } else {
        equalities["edges"] &= e.isEqual(e2);
      }
    }
    equalities["edgecount"] = m1.EdgeCount() == m2.EdgeCount();
    equalities["faces"] = true;
    for (face f = m1.faces(); f != m1.faces_end(); f++) {
      face f2 = m2.find_face(f());
      if (f2 == m2.faces_end()) {
        equalities["faces"] = false;
      } else {
        equalities["faces"] &= f.isEqual(f2);
      }
    }
    equalities["facecount"] = m1.FaceCount() == m2.FaceCount();


    /*for(face f=m1.faces(); f!=m1.faces_end(); f++){
      mout << f() << f.Left() << f.Right() << endl;
    }
    mout << "---------------------------" << endl;
    for(face f=m2.faces(); f!=m2.faces_end(); f++){
      mout << f() << f.Left() << f.Right() << endl;
    }*/
    bool res = true;
    for (auto &[key, equality] : equalities) {
      if (!equality) {
        mout << key << " is/are different" << endl;
        res = false;
      }
    }

    return res;
  }

  void TearDown() override { PPM->Barrier(0); }
};

INSTANTIATE_TEST_SUITE_P(TestRefineMesh, TestRefineMesh,
                         Values(TestConfig{"Interval2Refined", "Interval3Refined", "Stripes"}));


#endif // SPACETIME_TESTREFINEMESH_HPP
