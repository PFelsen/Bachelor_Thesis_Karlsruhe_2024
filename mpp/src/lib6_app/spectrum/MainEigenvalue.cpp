#include "m++.hpp"

#include "spectrum/EigenvalueExample.hpp"

/**
 * Example problem for eigenvalue methods on domain $\Omega:=(0,1)^2$
 *
 * $H:=H^1_0(\Omega)$
 * $M(u,v):=\int_\Omega \nabla u\cdot\nabla v - 4 x (1 - x) u v \,d(x,y)$
 * $N(u,v):=\int_\Omega u v \,d(x,y)$
 * -> consider eigenvalue problem: $M(u,v) = \lambda N(u,v)$ $(v\in H)$
 *
 * For the base problem we consider
 * $M_0(u,v):=\int_\Omega \nabla u\cdot\nabla v -  u v \,d(x,y)$
 * $N_0(u,v):=\int_\Omega u v \,d(x,y)$
 *
 * For the eigenvalue computations an additional shift parameter $\nu>0$ is introduced:
 * $M_\nu(u,v):=\int_\Omega \nabla u\cdot\nabla v + (\nu- 4 x (1 - x)) u v \,d(x,y)$
 * $N(u,v):=\int_\Omega u v \,d(x,y)$
 * -> shifted eigenvalue problem: $M_\nu(u,v) = (\lambda + \nu) N(u,v)$ $(v\in H)$
 *
 * Bilinear form considered during the homotpy
 * $M_t(u,v):=\int_\Omega \nabla u\cdot\nabla v + (\nu - 1 + t (1 - 4 x (1 - x))) u v \,d(x,y)$
 *
 * Goerisch's XbT-setting:
 * $X:=H(div, \Omega)\times L^2(\Omega)\times L^2(\Omega)$
 * $b_t(w,v):=\int_\Omega w_1\cdot v_1\,d(x,y) + \int_\Omega (\nu - 2 + t (1 - 4 x (1 - x))) w_2 v_2 \,d(x,y) + \int_\Omega w_3 v_3 \,d(x,y)$
 * $T:H\to X, Tu:=(\nabla u, u, u)$
 * Then $b_t(Tu,Tv)=M_t(u,v)$ holds true for all $u,v\in H$
 * Moreover, for $v\in H$ we can solve for the 3rd component:
 * -> $w_3:=v - (\nu - 2 + t (1 - 4 x (1 - x))) w_2 + div(w_1)$
 *
 * To compute $(w_1,w_2)\in H(div, \Omega)\times L^2(\Omega)$ corresponding to $v$ we approximately minimize
 * $b_t(w,w)=\int_\Omega |w_1|^2\,d(x,y) + \int_\Omega (\nu - 2 + t (1 - 4 x (1 - x))) |w_2|^2 \,d(x,y) + \int_\Omega |v - (\nu - 2 + t (1 - 4 x (1 - x))) w_2 + div(w_1)|^2 \,d(x,y)$
 */

int main(int argc, char **argv) {
  Config::SetConfPath(std::string(ProjectMppDir) + "/conf/");
  Mpp::initialize(&argc, argv);

  mout.setFileEnabled(false);

  // Running verified example
  RunExample<true>();

  return 0;
}
