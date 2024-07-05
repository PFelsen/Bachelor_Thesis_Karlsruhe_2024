#ifndef TUTORIAL_PLATEPROBLEM_HPP
#define TUTORIAL_PLATEPROBLEM_HPP

#include "IProblem.hpp"

class PlateProblem : public IProblem {
public:
  std::list<Point> corners = {Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0)};  
  PlateProblem () : IProblem() {
    std::string meshname = "UnitSquare4Triangles";
    Config::Get("Mesh", meshname);
    if (meshname == "Cook19Triangles")
      corners = {Point(0.0, 0.0), Point(48.0, 44.0), Point(48.0, 60.0), Point(0.0, 44.0)};  
    else if (meshname != "UnitSquare4Triangles")
      Exit("corners required for the boundary conditions");
  }
  string Name() const {
    return "PlateProblem";
  }
  double InitialValue (const Point& x) const {
    return 0;
  }
  double DirichletValue (const Point& x) const { return 0; }
  double Load (const Point& x) const { return 1; }
};

PlateProblem *CreatePlateProblem(const string &problemName);

std::shared_ptr<PlateProblem> CreatePlateProblemShared(const string &problemName);

#endif // TUTORIAL_PLATEPROBLEM_HPP
