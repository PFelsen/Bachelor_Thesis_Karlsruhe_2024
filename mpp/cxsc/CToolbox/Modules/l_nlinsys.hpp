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

/* CVS $Id: nlinsys.hpp,v 1.16 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: nlinsys (header)
// Purpose: Computing enclosures for all solutions of systems of nonlinear
//    equations given by continuously differentiable functions.
// Global functions:
//    AllNLSS()      : computes enclosures for all solutions
//    AllNLSSErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef _L_NLINSYS_HPP
#define _L_NLINSYS_HPP

#include <intvector.hpp>    // Integer vector type
#include <l_imatrix.hpp>      // Interval matrix/vector arithmetic
//#include <mv_util.hpp>      // Real matrix/vector utility functions
//#include <mvi_util.hpp>     // Interval matrix/vector utility functions
#include <l_grad_ari.hpp>     // Gradient differentiation arithmetic


using namespace cxsc;
using namespace std;


const int l_MaxCountNLSS = 10000; // Maximum count of result components

extern void NewtonStep( l_GTvector_FctPtr f, l_rvector y );

extern char* l_AllNLSSErrMsg ( int );
extern void  l_AllNLSS ( l_GTvector_FctPtr, l_ivector, real, l_imatrix&, intvector&,
                       int&, int&, int = l_MaxCountNLSS );
#endif




