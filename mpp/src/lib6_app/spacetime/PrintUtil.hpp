#ifndef SPACETIME_PRINTUTIL_H
#define SPACETIME_PRINTUTIL_H

#include "RMatrix.hpp"
#include "RVector.hpp"
#include <fstream>
#include "IMatrixGraph.hpp"
#include "Sparse.hpp"


std::string format_double(const double v);

std::string to_string(const std::vector<int> &V);

std::string to_string(const RVector &V);

void PrintScalar(std::ostream &s, const Scalar &a);

std::ostream &operator<<(std::ostream &s, const vector<int> &v);

void PrintRVectorHorizontal(std::ostream &s, const RVector &v);

std::ostream &operator<<(std::ostream &s, const std::pair<RVector, double> &va);

void PrintValues2(const std::string &name, const RVector &V);

void PrintValues(const std::string &name, const RVector &V, int run = 100, int verbose = -1);

void PrintMemoryInfo(const IMatrixGraph &mg_fine);

#endif //SPACETIME_PRINTUTIL_H
