#ifndef SPACETIMEMEASUREMENT_HPP
#define SPACETIMEMEASUREMENT_HPP

#include "SpaceTime.hpp"
#include "STAssemble.hpp"
#include "AcousticProblems.hpp"
#include <memory>

void measurePkt(const Vector &u,
                const AcousticProblem &problem,
                const STDiscretization &disc,
                int run);


#endif //SPACETIMEMEASUREMENT_HPP
