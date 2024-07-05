#ifndef TESTEIGENSOLVER_HPP
#define TESTEIGENSOLVER_HPP

#include "EigenSolverCreator.hpp"
#include "MeshesCreator.hpp"
#include "TestEnvironment.hpp"

struct EigenSolverTestParameter {
  std::string meshName;
  int level;
  string discName;
  int degree;
  std::list<Point> corners; // only needed for Argyris
};

static vector<EigenSolverTestParameter> TestCases = {
    EigenSolverTestParameter{"Square2Triangles", 6, "Lagrange", 2, {}},
    EigenSolverTestParameter{"Square2Triangles", 5, "Lagrange", 3, {}},
    EigenSolverTestParameter{"Square", 6, "Serendipity", 2, {}},
    EigenSolverTestParameter{"Square", 5, "Lagrange", 3, {}},
    EigenSolverTestParameter{"SquareTurned", 5, "Lagrange", 3, {}},
    EigenSolverTestParameter{"Square4Triangles",
                             5,
                             "Argyris",
                             -1,
                             {Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 1.0),
                              Point(0.0, 1.0)}} //,
    //        EigenSolverTestParameter{
    //                "Square4TrianglesTurned", 6, "Argyris", -1,
    //                {Point(0.70710678118, 0.0), Point(0.0,0.70710678118), Point(-0.70710678118,
    //                0.0), Point(0.0, -0.70710678118)}
    //        }
};


#endif // TESTEIGENSOLVER_HPP
