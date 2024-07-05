#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include "QuadratureGL.hpp"
#include "QuadratureSym.hpp"

template<typename T, int sDim, int tDim>
const QuadratureT<T, sDim, tDim> &GetQuadratureT(const CELLTYPE, int exactUpTo = -1,
                                                 int exactUpTo1 = -1, int exactUpTo2 = -1);

const Quadrature &GetQuadrature(const CELLTYPE, int exactUpTo = -1, int exactUpTo1 = -1,
                                int exactUpTo2 = -1);

template<typename T, int sDim, int tDim>
const QuadratureT<T, sDim, tDim> &GetGLQuadratureT(const CELLTYPE, int npcount);

const Quadrature &GetGLQuadrature(const CELLTYPE, int npcount);

template<class T, int sDim, int tDim>
const QuadratureT<T, sDim, tDim> &GetQuadrature(const std::string &name);

const Quadrature &GetQuadrature(const std::string &);

extern void clearQuadrature();

#endif // of #ifndef _QUADRATURE_H_
