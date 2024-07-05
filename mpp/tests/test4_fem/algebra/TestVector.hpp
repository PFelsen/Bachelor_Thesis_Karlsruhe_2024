#ifndef FWI_TESTVECTOR_HPP
#define FWI_TESTVECTOR_HPP

#include <memory>

#include <LagrangeDiscretization.hpp>
#include <Meshes.hpp>
#include <MeshesCreator.hpp>
#include <Vector.hpp>

#include "LevelPair.hpp"
#include "TestEnvironment.hpp"

class TestVector : public Test {
protected:
  std::string meshName = "Interval";

  std::string distName = "RCB";

  int commSplit = 0;

  int pLevel = 2;

  int level = 3;

  std::shared_ptr<const Meshes> meshes;
  std::shared_ptr<const LagrangeDiscretization> disc;
  Vector vector;
public:
  TestVector() :
      meshes(MeshesCreator(meshName)
                 .WithDistribute(distName)
                 .WithPLevel(pLevel)
                 .WithLevel(level)
                 .Create()),
      disc(std::make_shared<const LagrangeDiscretization>(*meshes, 0, 1)),
      vector(Vector(1.0, disc)) {
    // Ensure meshes are generated
    (*meshes)[LevelPair{level, -1, 0, commSplit}];

    meshes->PrintInfo();
  }

  void testVectorNonOverlap(const Vector &vector, double testval) {
    std::vector<bool> equal;
    for (row r = vector.rows(); r != vector.rows_end(); ++r) {
      procset p = vector.find_procset(r());
      if (p == vector.procsets_end()) {
        equal.push_back(vector(r, 0) == testval);
      } else {
        continue;
      }
    }
    std::vector<bool> equaltrue(equal.size(), true);
    EXPECT_VECTOR_EQ(equaltrue, equal);
  }

  void testVectorOverlap(const Vector &vector, double testval) {
    std::vector<bool> equal;
    for (row r = vector.rows(); r != vector.rows_end(); ++r) {
      procset p = vector.find_procset(r());
      if (p == vector.procsets_end()) continue;
      else {
        size_t psize = vector.find_procset(r()).size();
        equal.push_back(vector(r, 0) == testval * psize);
      }
    }
    std::vector<bool> equaltrue(equal.size(), true);
    EXPECT_VECTOR_EQ(equaltrue, equal);
  }

  bool addRowWithSingleOwner(Vector &v, const Vector &u, const procset &p, const row &r) {
    if (p == u.procsets_end()) {
      for (int i = 0; i < r.n(); ++i) {
        v(r, i) += u(r, i);
      }
      return true;
    } else return false;
  };

  void addOnlyOnMaster(Vector &v, const Vector &u, const procset &p, const row &r) {
    if (p.master() == PPM->Proc(0)) {
      for (int i = 0; i < r.n(); ++i) {
        v(r, i) += u(r, i);
      }
    }
  };

  bool addNormFromSingleOwner(const Vector &v, const procset &p, const row &r, double &norm) {
    if (p == v.procsets_end()) {
      for (int i = 0; i < r.n(); ++i) {
        norm += v(r, i) * v(r, i);
      }
      return true;
    } else return false;
  };

  void addNormOnMaster(const Vector &v, const procset &p, const row &r, double &norm) {
    if (p.master() == PPM->Proc(0)) {
      for (int i = 0; i < r.n(); ++i) {
        norm += v(r, i) * v(r, i);
      }
    }
  };

  void testAdditionWithoutOverlap(Vector &v, const Vector &u) {
    for (row r = u.rows(); r != u.rows_end(); ++r) {
      procset p = u.find_procset(r());
      if (addRowWithSingleOwner(v, u, p, r)) continue;
      else addOnlyOnMaster(v, u, p, r);
    }
  }

  void testNormWithoutOverlap(Vector &v) {
    double norm = 0.0;
    for (row r = v.rows(); r != v.rows_end(); ++r) {
      procset p = v.find_procset(r());
      if (addNormFromSingleOwner(v, p, r, norm)) continue;
      else addNormOnMaster(v, p, r, norm);
    }
    norm = sqrt(PPM->SumOnCommSplit(norm, 0));
    EXPECT_EQ(norm, sqrt(pow(2, level)));
  }

  void TearDown() override { PPM->Barrier(0); }
};

#endif // FWI_TESTVECTOR_HPP
