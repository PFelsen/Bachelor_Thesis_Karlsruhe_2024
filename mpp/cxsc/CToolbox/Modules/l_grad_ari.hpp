/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: grad_ari.hpp,v 1.16 2014/01/30 17:49:26 cxsc Exp $ */

//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz   
//
// For details on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
//============================================================================
//----------------------------------------------------------------------------
//
//       Modified for C-XSC++ (Version 0.9.1) Stefan Wedner May 2000
//
//----------------------------------------------------------------------------
//
// File: grad_ari (header)
// Purpose: Definition of a multi-dimensional l_interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to first order (i.e. gradient or Jacobian matrix).
// Types:
//    l_GTscalar_FctPtr         : pointer for a scalar valued function
//    l_GTvector_FctPtr         : pointer for a vector valued function
// Class l_GradType:
//    l_GradType()              : constructors
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for arguments of types
//                              l_GradType, l_interval, and real
//    GradVar()               : to define l_GradType variables
//    fValue(), gradValue()   : to get function and derivative values
//    operators +, -, *, /    : operators of diff. arithmetic
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEvalG()                : to compute function value only
//    fgEvalG()               : to compute function and first derivative
//                              value (gradient)
// Class l_GTvector:
//    l_GTvector()              : constructors
//    ~l_GTvector()             : destructor
//    Dim()                   : to get the actual dimension
//    operator =              : assignment operator
//    operator []             : component access
//    fValue(), JacValue()    : to get function and derivative values
//    fEvalJ()                : to compute function value only
//    fJEvalJ()               : to compute function and first derivative
//                              value (Jacobian matrix)
//----------------------------------------------------------------------------
#ifndef _L_GRAD_ARI_HPP
#define _L_GRAD_ARI_HPP

#include <l_interval.hpp>     // l_interval arithmetic
#include <l_imatrix.hpp>      // l_interval matrix/vector arithmetic
#include <mvi_util.hpp>     // l_interval matrix/vector utility functions

using namespace cxsc;
using namespace std;

class l_GradType;
class l_GTvector;

typedef l_GradType (*l_GTscalar_FctPtr)(const l_GTvector&);
typedef l_GTvector (*l_GTvector_FctPtr)(const l_GTvector&);

class l_GradType {
    int     nmax;
    l_ivector g;

    friend void TestSize ( const l_GradType&, const l_GradType&, const char* );

  public:
    l_GradType ( );
    explicit l_GradType ( int );
    l_GradType ( const l_GradType& );

    friend void Resize ( l_GradType&, int );
    int  Dim ( ) const { return nmax; };

    l_interval& operator[] ( int );

    l_GradType& operator= ( const l_GradType& );
    l_GradType& operator= ( const l_interval& );
    l_GradType& operator= ( const real& );


    friend l_GTvector GradVar ( const l_ivector& );
    friend l_GTvector GradVar ( const l_rvector& );

    friend l_interval fValue    ( const l_GradType& );
    friend l_ivector  gradValue ( const l_GradType& );

    friend l_GradType operator+ ( l_GradType&);
    friend l_GradType operator- ( const l_GradType&);
    friend l_GradType operator+ ( const l_GradType&, const l_GradType&);
    friend l_GradType operator- ( const l_GradType&, const l_GradType&);
    friend l_GradType operator* ( const l_GradType&, const l_GradType&);
    friend l_GradType operator/ ( const l_GradType&, const l_GradType&);
    friend l_GradType operator+ ( const l_GradType&, const l_interval&);
    friend l_GradType operator- ( const l_GradType&, const l_interval&);
    friend l_GradType operator* ( const l_GradType&, const l_interval&);
    friend l_GradType operator/ ( const l_GradType&, const l_interval&);
    friend l_GradType operator+ ( const l_interval&, const l_GradType&);
    friend l_GradType operator- ( const l_interval&, const l_GradType&);
    friend l_GradType operator* ( const l_interval&, const l_GradType&);
    friend l_GradType operator/ ( const l_interval&, const l_GradType&);
    friend l_GradType operator+ ( const l_GradType&, const real&);
    friend l_GradType operator- ( const l_GradType&, const real&);
    friend l_GradType operator* ( const l_GradType&, const real&);
    friend l_GradType operator/ ( const l_GradType&, const real&);
    friend l_GradType operator+ ( const real&, const l_GradType&);
    friend l_GradType operator- ( const real&, const l_GradType&);
    friend l_GradType operator* ( const real&, const l_GradType&);
    friend l_GradType operator/ ( const real&, const l_GradType&);
    friend l_GradType sqr   ( const l_GradType& );
    friend l_GradType power ( const l_GradType&, const int );
    friend l_GradType sqrt  ( const l_GradType& );
    friend l_GradType exp   ( const l_GradType& );
    friend l_GradType ln    ( const l_GradType& );
    friend l_GradType sin   ( const l_GradType& );
    friend l_GradType cos   ( const l_GradType& );
    friend l_GradType tan   ( const l_GradType& );
    friend l_GradType cot   ( const l_GradType& );
    friend l_GradType asin  ( const l_GradType& );
    friend l_GradType acos  ( const l_GradType& );
    friend l_GradType atan  ( const l_GradType& );
    friend l_GradType acot  ( const l_GradType& );
    friend l_GradType sinh  ( const l_GradType& );
    friend l_GradType cosh  ( const l_GradType& );
    friend l_GradType tanh  ( const l_GradType& );
    friend l_GradType coth  ( const l_GradType& );
    friend l_GradType asinh ( const l_GradType& );
    friend l_GradType acosh ( const l_GradType& );
    friend l_GradType atanh ( const l_GradType& );
    friend l_GradType acoth ( const l_GradType& );

    friend void fEvalG   ( l_GTscalar_FctPtr, l_ivector, l_interval& );
    friend void fgEvalG  ( l_GTscalar_FctPtr, l_ivector, l_interval&, l_ivector& );
};

//============================================================================

void TestSize ( const l_GradType&, const l_GradType&, const char* );
void Resize ( l_GradType&, int );
l_GTvector GradVar ( const l_ivector& );
l_GTvector GradVar ( const l_rvector& );

l_interval fValue    ( const l_GradType& );
l_ivector  gradValue ( const l_GradType& );

l_GradType operator+ ( l_GradType&);
l_GradType operator- ( const l_GradType&);
l_GradType operator+ ( const l_GradType&, const l_GradType&);
l_GradType operator- ( const l_GradType&, const l_GradType&);
l_GradType operator* ( const l_GradType&, const l_GradType&);
l_GradType operator/ ( const l_GradType&, const l_GradType&);
l_GradType operator+ ( const l_GradType&, const l_interval&);
l_GradType operator- ( const l_GradType&, const l_interval&);
l_GradType operator* ( const l_GradType&, const l_interval&);
l_GradType operator/ ( const l_GradType&, const l_interval&);
l_GradType operator+ ( const l_interval&, const l_GradType&);
l_GradType operator- ( const l_interval&, const l_GradType&);
l_GradType operator* ( const l_interval&, const l_GradType&);
l_GradType operator/ ( const l_interval&, const l_GradType&);
l_GradType operator+ ( const l_GradType&, const real&);
l_GradType operator- ( const l_GradType&, const real&);
l_GradType operator* ( const l_GradType&, const real&);
l_GradType operator/ ( const l_GradType&, const real&);
l_GradType operator+ ( const real&, const l_GradType&);
l_GradType operator- ( const real&, const l_GradType&);
l_GradType operator* ( const real&, const l_GradType&);
l_GradType operator/ ( const real&, const l_GradType&);
l_GradType sqr   ( const l_GradType& );
l_GradType power ( const l_GradType&, const int );
l_GradType sqrt  ( const l_GradType& );
l_GradType exp   ( const l_GradType& );
l_GradType ln    ( const l_GradType& );
l_GradType sin   ( const l_GradType& );
l_GradType cos   ( const l_GradType& );
l_GradType tan   ( const l_GradType& );
l_GradType cot   ( const l_GradType& );
l_GradType asin  ( const l_GradType& );
l_GradType acos  ( const l_GradType& );
l_GradType atan  ( const l_GradType& );
l_GradType acot  ( const l_GradType& );
l_GradType sinh  ( const l_GradType& );
l_GradType cosh  ( const l_GradType& );
l_GradType tanh  ( const l_GradType& );
l_GradType coth  ( const l_GradType& );
l_GradType asinh ( const l_GradType& );
l_GradType acosh ( const l_GradType& );
l_GradType atanh ( const l_GradType& );
l_GradType acoth ( const l_GradType& );

void fEvalG   ( l_GTscalar_FctPtr, l_ivector, l_interval& );
void fgEvalG  ( l_GTscalar_FctPtr, l_ivector, l_interval&, l_ivector& );

//============================================================================

class l_GTvector {
    int      nmax;
    l_GradType *gt;

    friend void TestSize ( const l_GTvector&, const l_GTvector&, const char* );

  public:
    explicit l_GTvector  ( int );
    l_GTvector  ( const l_GTvector& );
    ~l_GTvector ( );

    int Dim ( ) const { return nmax; };

    l_GTvector& operator=  ( const l_GTvector& );

    l_GradType& operator[] ( int );
    const l_GradType& operator[] ( int ) const;

    friend l_ivector  fValue   ( const l_GTvector& );
    friend l_imatrix  JacValue ( const l_GTvector& );

    friend void fEvalJ  ( l_GTvector_FctPtr, l_ivector, l_ivector& );
    friend void fJEvalJ ( l_GTvector_FctPtr, l_ivector, l_ivector&, l_imatrix& );
};

//============================================================================

void TestSize ( const l_GTvector&, const l_GTvector&, const char* );

l_ivector  fValue   ( const l_GTvector& );
l_imatrix  JacValue ( const l_GTvector& );

void fEvalJ  ( l_GTvector_FctPtr, l_ivector, l_ivector& );
void fJEvalJ ( l_GTvector_FctPtr, l_ivector, l_ivector&, l_imatrix& );

#endif




