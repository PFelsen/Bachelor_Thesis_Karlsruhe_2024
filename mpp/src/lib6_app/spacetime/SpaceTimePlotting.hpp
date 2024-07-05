#ifndef SPACETIMEPLOTTING_HPP
#define SPACETIMEPLOTTING_HPP

#include "SpaceTimeTools.hpp"
#include "STAssemble.hpp"
#include "Shapes.hpp"

void printSingleSolutionVTK(const Vector &U,
                            const STAssemble &assemble,
                            const string &name);


void printVTK(const Vector &U,
              const STAssemble &assemble,
              int run = 0);

void printVTK_Eta(const Vector &Eta,
                  const STAssemble &assemble,
                  const std::string &filename);

void printVTK_dual(const Vector &U_dual,
                   const Vector &Eta,
                   const STAssemble &assemble,
                   int run);


#endif //SPACETIMEPLOTTING_HPP
