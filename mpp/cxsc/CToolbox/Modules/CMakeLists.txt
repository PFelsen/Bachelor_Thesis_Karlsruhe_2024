set(ctoolbox_src
  ci_util.cpp
  cipoly.cpp
  cpoly.cpp
  cpzero.cpp
  ddf_ari.cpp
  expreval.cpp
  fastlss.cpp
  gop.cpp
  gop1.cpp
  grad_ari.cpp
  hess_ari.cpp
  linsys.cpp
  lop.cpp
  lop_ari.cpp
  lst1_ari.cpp
  lst_ari.cpp
  matinv.cpp
  nlfzero.cpp
  nlinsys.cpp
  rev_simp.cpp
  rpeval.cpp
  rpoly.cpp
  set_ari.cpp
  xi_ari.cpp
  l_grad_ari.cpp
  l_matinv.cpp
  l_xi_ari.cpp
  l_nlinsys.cpp
)

include_directories(${PROJECT_MPP_DIR}/cxsc/CToolbox/Modules)

add_library(cxsc_toolbox STATIC ${ctoolbox_src})

target_link_libraries(cxsc_toolbox cxsc)