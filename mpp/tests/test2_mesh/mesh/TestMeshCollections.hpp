#ifndef TESTMESHCOLLECTIONS_H
#define TESTMESHCOLLECTIONS_H


#include "TestEnvironment.hpp"
#include "TestMesh.hpp"

#include "CoarseGeometry.hpp"
#include "Mesh.hpp"


static MeshTestParameters testMeshParams2D(std::vector<vector<double>>{{0.0, 0.0, 0.0},  // 0
                                                                       {1.0, 0.0, 0.0},  // 1
                                                                       {2.0, 0.0, 0.0},  // 2
                                                                       {3.0, 0.0, 0.0},  // 3
                                                                       {4.0, 0.0, 0.0},  // 4
                                                                       {4.0, 1.0, 0.0},  // 5
                                                                       {3.0, 1.0, 0.0},  // 6
                                                                       {2.0, 1.0, 0.0},  // 7
                                                                       {1.0, 1.0, 0.0},  // 8
                                                                       {0.0, 1.0, 0.0}}, // 9
                                           std::vector<vector<int>>{{4, 10, 0, 1, 8, 9},
                                                                    {4, 10, 1, 2, 7, 8},
                                                                    {4, 10, 2, 3, 6, 7},
                                                                    {4, 10, 3, 4, 5, 6}},
                                           std::vector<vector<int>>{{2, 5, 0, 1},
                                                                    {2, 5, 1, 2},
                                                                    {2, 5, 2, 3},
                                                                    {2, 5, 3, 4},
                                                                    {2, 5, 4, 5},
                                                                    {2, 5, 5, 6},
                                                                    {2, 5, 6, 7},
                                                                    {2, 5, 7, 8},
                                                                    {2, 5, 8, 9},
                                                                    {2, 5, 9, 0}});
static vector<MeshTestParameters>
    unitGeometries{// Unit Interval
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                      std::vector<vector<int>>{{2, 10, 0, 1}},
                                      std::vector<vector<int>>{{1, 5, 0}, {1, 15, 1}}}
#if SpaceDimension >= 2
                   ,
                   // Unit Triangle
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {0.0, 1.0, 0.0}},
                                      std::vector<vector<int>>{{3, 10, 0, 1, 2}},
                                      std::vector<vector<int>>{{2, 5, 0, 1},
                                                               {2, 15, 1, 2},
                                                               {2, 25, 2, 0}}},
                   // Unit Square
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {1.0, 1.0, 0.0},
                                                                  {0.0, 1.0, 0.0}},
                                      std::vector<vector<int>>{{4, 10, 0, 1, 2, 3}},
                                      std::vector<vector<int>>{{2, 5, 0, 1},
                                                               {2, 15, 1, 2},
                                                               {2, 25, 2, 3},
                                                               {2, 35, 3, 0}}}
#endif
#if SpaceDimension >= 3
                   ,
                   // Unit Tetrahedron
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {0.0, 1.0, 0.0},
                                                                  {0.0, 0.0, 1.0}},
                                      std::vector<vector<int>>{{4, 10, 0, 1, 2, 3}},
                                      std::vector<vector<int>>{{3, 5, 0, 1, 2},
                                                               {3, 15, 0, 1, 3},
                                                               {3, 25, 1, 2, 3},
                                                               {3, 35, 2, 0, 3}}},
                   // Unit Hexahedron
                   MeshTestParameters{std::vector<vector<double>>{{0.0, 0.0, 0.0},
                                                                  {1.0, 0.0, 0.0},
                                                                  {1.0, 1.0, 0.0},
                                                                  {0.0, 1.0, 0.0},
                                                                  {0.0, 0.0, 1.0},
                                                                  {1.0, 0.0, 1.0},
                                                                  {1.0, 1.0, 1.0},
                                                                  {0.0, 1.0, 1.0}},
                                      std::vector<vector<int>>{{8, 10, 0, 1, 2, 3, 4, 5, 6, 7}},
                                      std::vector<vector<int>>{{4, 5, 0, 1, 2, 3},
                                                               {4, 15, 0, 1, 5, 4},
                                                               {4, 25, 1, 2, 6, 5},
                                                               {4, 35, 2, 3, 7, 6},
                                                               {4, 45, 3, 0, 4, 7},
                                                               {4, 55, 4, 5, 6, 7}}}
#endif
    };

#endif // TESTMESHCOLLECTIONS_H
