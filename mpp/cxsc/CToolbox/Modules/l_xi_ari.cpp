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

/* CVS $Id: xi_ari.cpp,v 1.16 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: xi_ari (implementation)
// Purpose: Definition of an extended l_interval arithmetic which allows the
//    division by an l_interval containing zero.
// Method: Overloading of operators for arithmetic and lattice operations
//    of data type 'l_xinterval'.
// Global type:
//    KindType         : component type of extended l_intervals
// Global function:
//    EmptyIntval()    : empty set as irregular l_interval
// Class l_xinterval:
//    operators %, -, &: operators of extended l_interval arithmetic
//    operator =       : assignment operator
// Updates:
//    04.03.1996 : Modification of extended l_interval operations and kind type
//----------------------------------------------------------------------------
//#include <i_util.hpp>     // l_interval utility functions
#include <l_xi_ari.hpp>

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// An extended l_interval 'x', represented by the type 'l_xinterval', is defined
// according to the following rules (a <= b):
//
// x = [a, b]              :  x.kind = Finite,     x.inf = a, x.sup = b
// x = [a, +oo]            :  x.kind = PlusInfty,  x.inf = a, x.sup undefined
// x = [-oo, a]            :  x.kind = MinusInfty, x.sup = a, x.inf undefined
// x = [-oo, a] v [b, +oo] :  x.kind = Double,     x.inf = b, x.sup = a
// x = [-oo, +oo]          :  x.kind = Double,     x.inf = a, x.sup = a
// x = [/]                 :  x.kind = Empty,      x.inf  and x.sup undefined
//
// In this definition, 'v' stands for the set union and 'oo' for infinity.
//----------------------------------------------------------------------------
l_interval l_EmptyIntval ( )                         // Irregular (empty) l_interval
{                                                //---------------------------
  l_interval x;
  Inf(x) = l_real(999999999.0);
  Sup(x) = l_real(-999999999.0);
  return x;
}

//----------------------------------------------------------------------------
// Constructors and assignment operator
//----------------------------------------------------------------------------
l_xinterval::l_xinterval ( )
{
  kind = Finite;
  inf  = 0.0;
  sup  = 0.0;
}

l_xinterval::l_xinterval ( const KindType& k, const l_real& i, const l_real& s )
{
  kind = k;
  inf  = i;
  sup  = s;
}

l_xinterval::l_xinterval ( const l_xinterval& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
}

l_xinterval& l_xinterval::operator= ( const l_xinterval& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
  return *this;
}

//----------------------------------------------------------------------------
// Extended l_interval division 'A / B' where 0 in 'B' is allowed.
//----------------------------------------------------------------------------
l_xinterval operator% ( const l_interval& A, const l_interval& B )
{
  l_interval  c;
  l_xinterval Q;

  if ( in(0.0, B) ) {
    if ( in(0.0, A) ) {
      Q.kind = Double;                    // Q = [-oo,+oo] = [-oo,0] v [0,+oo]
      Q.sup  = 0.0;                       //----------------------------------
      Q.inf  = 0.0;
    }
    else if ( B == 0.0 ) {                                          // Q = [/]
      Q.kind = PlusInfty;                                           //--------
      Q.inf  = Inf(A/B);//  divd(Sup(A),Inf(B));
    }
    else if ( (Sup(A) < 0.0) && (Sup(B) == 0.0) ) {         // Q = [Q.inf,+oo]
      Q.kind = PlusInfty;                                   //----------------
      Q.inf  = Inf(A/B);//  divd(Sup(A),Inf(B));
    }
    else if ( (Sup(A) < 0.0) && (Inf(B) < 0.0) && (Sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      Q.sup  = Sup(A/B);//  divu(Sup(A),Sup(B));           //------------------------------
      Q.inf  = Inf(A/B);//  divd(Sup(A),Inf(B));
    }
    else if ( (Sup(A) < 0.0) && (Inf(B) == 0.0) ) {         // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      Q.sup  = Sup(A/B);//  divu(Sup(A),Sup(B));
    }
    else if ( (Inf(A) > 0.0) && (Sup(B) == 0.0) ) {         // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      Q.sup  = Sup(A/B);//  divu(Inf(A),Inf(B));
    }
    else if ( (Inf(A) > 0.0) && (Inf(B) < 0.0) && (Sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      Q.sup  = Sup(A/B);//  divu(Inf(A),Inf(B));           //------------------------------
      Q.inf  = Inf(A/B);//  divd(Inf(A),Sup(B));
    }
    else { // if ( (Inf(A) > 0.0) && (Inf(B) == 0.0) )
      Q.kind = PlusInfty;                                   // Q = [Q.inf,+oo]
      Q.inf  = Inf(A/B);//  divd(Inf(A),Sup(B));                         //----------------
    }
  } // in(0.0,B)
  else {  // !in(0.0,B)
    c = A / B;                                            // Q = [C.inf,C.sup]
    Q.kind = Finite;                                      //------------------
    Q.inf  = Inf(c);
    Q.sup  = Sup(c);
  }

  return Q;
} // operator%

//----------------------------------------------------------------------------
// Subtraction of an extended l_interval 'B' from a l_real value 'a'.
//----------------------------------------------------------------------------
l_xinterval operator- ( const l_real& a, const l_xinterval& B )
{
  l_xinterval D;

  switch (B.kind) {
    case Finite     : D.kind = Finite;                    // D = [D.inf,D.sup]
                      D.inf  = Inf(a-l_interval(B.sup));// subd(a,B.sup);             //------------------
                      D.sup  = Sup(a-l_interval(B.inf));// subu(a,B.inf);
                      break;
    case PlusInfty  : D.kind = MinusInfty;                    // D = [inf,+oo]
                      D.sup  = Sup(a-l_interval(B.inf));// subu(a,B.inf);                 //--------------
                      break;
    case MinusInfty : D.kind = PlusInfty;                     // D = [-oo,sup]
                      D.inf  = Inf(a-l_interval(B.sup));// subd(a,B.sup);                 //--------------
                      break;
    case Double     : D.kind = Double;        // D = [-oo,D.sup] v [D.inf,+oo]
                      D.inf  = Inf(a-l_interval(B.sup));// subd(a,B.sup); //------------------------------
                      D.sup  = Sup(a-l_interval(B.inf));// subu(a,B.inf);
                      if (D.inf < D.sup) D.inf = D.sup;
                      break;
    case Empty      : D.kind = Empty;                               // D = [/]
                      D.inf  = Inf(a-l_interval(B.sup));// subd(a,B.sup);                       //--------
                      break;
  } // switch
  return D;
}

//----------------------------------------------------------------------------
// Intersection of an l_interval 'X' and an extended l_interval 'Y'. The result
// is given as a pair (vector) of l_intervals, where one or both of them can
// be empty l_intervals.
//----------------------------------------------------------------------------
l_ivector operator& ( const l_interval& X, const l_xinterval& Y )
{
  l_interval H;
  l_ivector  IS(2);

  IS[1] = EmptyIntval();
  IS[2] = EmptyIntval();

  switch (Y.kind) {
    case Finite     : // [X.inf,X.sup] & [Y.inf,Y.sup]
                      //------------------------------
                      H = l_interval(Y.inf,Y.sup);
                      if ( !Disjoint(X,H) ) IS[1] = X & H;
                      break;
    case PlusInfty  : // [X.inf,X.sup] & [Y.inf,+oo]
                      //----------------------------
                      if (Sup(X) >= Y.inf) {
                        if (Inf(X) > Y.inf)
                          IS[1] = X;
                        else
                          IS[1] = l_interval(Y.inf,Sup(X));
		      }  
                      break;
    case MinusInfty : // [X.inf,X.sup] & [-oo,Y.sup]
                      //----------------------------
                      if (Y.sup >= Inf(X)) {
                        if (Sup(X)<Y.sup)
                          IS[1] = X;
                        else
                          IS[1] = l_interval(Inf(X),Y.sup);
		      }	  
                      break;
    case Double     : if ( (Inf(X) <= Y.sup) && (Y.inf <= Sup(X)) ) {
                        IS[1] = l_interval(Inf(X),Y.sup);    // X & [-oo,Y.sup]
                        IS[2] = l_interval(Y.inf,Sup(X));    // X & [Y.inf,+oo]
                      }
                      else if (Y.inf <= Sup(X)) // [X.inf,X.sup] & [Y.inf,+oo]
                        if (Inf(X) >= Y.inf)    //----------------------------
                          IS[1] = X;
                        else
                          IS[1] = l_interval(Y.inf,Sup(X));
                      else if (Inf(X) <= Y.sup) { // [X.inf,X.sup] & [-oo,Y.sup]
                        if (Sup(X) <= Y.sup)      //----------------------------
                          IS[1] = X;
                        else
                          IS[1] = l_interval(Inf(X),Y.sup);
		      }
                      break;
    case Empty      : break;                           // [X.inf,X.sup] ** [/]
  } // switch                                          //---------------------

  return IS;
} // operator&




