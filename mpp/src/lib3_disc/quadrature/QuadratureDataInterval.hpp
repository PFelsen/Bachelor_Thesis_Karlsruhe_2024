#ifndef QUADRATUREDATAINTERVAL_HPP
#define QUADRATUREDATAINTERVAL_HPP

#ifdef BUILD_IA
#include "IAInterval.hpp"
#endif

/**
 * The verification files for verified quadrature rules using interval arithmetic are given in
 * test3_disc/quadrature/verifyquadrature!
 */

template<typename T = double>
struct QintSymT_data0 {
  T weight;
};

template<typename T = double>
struct QintSymT_data1 {
  T weight;
  T A;
  T B;
};

template<typename T = double>
struct QintSymT_data {
  std::vector<QintSymT_data0<T>> d0;
  std::vector<QintSymT_data1<T>> d1;

  int Size() { return d0.size() + 2 * d1.size(); }
};

template<typename T = double>
QintSymT_data<T> GetQintSymT_data(int exactUpTo);

template<>
QintSymT_data<double> GetQintSymT_data<double>(int exactUpTo) {
  switch (exactUpTo) {
  case 0:
  case 1:
    return QintSymT_data<double>({{1.0}}, {});
  case 2:
  case 3:
    return QintSymT_data<double>({},
                                 {{0.5,
                                   0.2113248654051871177454256097490212721761991243649365619906,
                                   0.7886751345948128822545743902509787278238008756350634380093}});
  case 4:
  case 5:
    return QintSymT_data<double>({{0.444444444444444}},
                                 {{0.277777777777777, 0.112701665379258, 0.887298334620742}});
  case 6:
  case 7:
    return QintSymT_data<double>({}, {{0.173927422568727, 0.069431844202974, 0.930568155797026},
                                      {0.326072577431273, 0.330009478207571, 0.669990521792429}});
  case 8:
  case 9:
    return QintSymT_data<double>({{0.284444444444444}},
                                 {{0.118463442528095, 0.046910077030667, 0.953089922969333},
                                  {0.239314335249683, 0.230765344947160, 0.769234655052840}});
  case 10:
  case 11:
    return QintSymT_data<double>({}, {{0.085662246189585, 0.033765242898424, 0.966234757101576},
                                      {0.180380786524069, 0.169395306766868, 0.830604693233132},
                                      {0.233956967286345, 0.380690406958402, 0.619309593041598}});
  case 12:
  case 13:
    return QintSymT_data<double>({{0.208979591836735}},
                                 {{0.064742483084435, 0.025446043828621, 0.974553956171379},
                                  {0.139852695744638, 0.129234407200303, 0.870765592799697},
                                  {0.190915025252560, 0.297077424311301, 0.702922575688699}});
  case 14:
  case 15:
    return QintSymT_data<double>({}, {{0.050614268145188, 0.019855071751232, 0.980144928248768},
                                      {0.111190517226687, 0.101666761293187, 0.898333238706813},
                                      {0.156853322938944, 0.237233795041836, 0.762766204958164},
                                      {0.181341891689181, 0.408282678752175, 0.591717321247825}});
  case 16:
  case 17:
    return QintSymT_data<double>({{0.165119677500630}},
                                 {{0.040637194180787, 0.015919880246187, 0.984080119753813},
                                  {0.090324080347429, 0.081984446336682, 0.918015553663318},
                                  {0.130305348201468, 0.193314283649705, 0.806685716350295},
                                  {0.156173538520001, 0.337873288298096, 0.662126711701905}});
  case 18:
  case 19:
    return QintSymT_data<double>({}, {{0.0333356721543440, 0.0130467357414141, 0.986953264258586},
                                      {0.0747256745752903, 0.0674683166555077, 0.932531683344492},
                                      {0.109543181257991, 0.160295215850488, 0.839704784149512},
                                      {0.134633359654998, 0.283302302935376, 0.716697697064624},
                                      {0.147762112357376, 0.425562830509184, 0.574437169490816}});
  default:
    Exit("Quadrature exact up to degree " + std::to_string(exactUpTo)
         + " not implemented for intervals!");
  }
}

#ifdef BUILD_IA

template<>
QintSymT_data<IAInterval> GetQintSymT_data<IAInterval>(int exactUpTo) {
  switch (exactUpTo) {
  case 0:
  case 1:
    return QintSymT_data<IAInterval>({{IAInterval(1.000000000000000, 1.000000000000000)}}, {});
  case 2:
  case 3:
    return QintSymT_data<IAInterval>({}, {{IAInterval(0.500000000000000, 0.500000000000000),
                                           IAInterval(0.211324865405187, 0.211324865405188),
                                           IAInterval(0.788675134594812, 0.788675134594813)}});
  case 4:
  case 5:
    return QintSymT_data<IAInterval>({{IAInterval(0.444444444444444, 0.444444444444445)}},
                                     {{IAInterval(0.277777777777777, 0.277777777777778),
                                       IAInterval(0.112701665379258, 0.112701665379259),
                                       IAInterval(0.887298334620741, 0.887298334620742)}});
  case 6:
  case 7:
    return QintSymT_data<IAInterval>({}, {{IAInterval(0.173927422568726, 0.173927422568727),
                                           IAInterval(0.069431844202973, 0.069431844202974),
                                           IAInterval(0.930568155797026, 0.930568155797027)},
                                          {IAInterval(0.326072577431273, 0.326072577431274),
                                           IAInterval(0.330009478207571, 0.330009478207572),
                                           IAInterval(0.669990521792428, 0.669990521792429)}});
  case 8:
  case 9:
    return QintSymT_data<IAInterval>({{IAInterval(0.284444444444444, 0.284444444444445)}},
                                     {{IAInterval(0.118463442528094, 0.118463442528095),
                                       IAInterval(0.046910077030667, 0.046910077030669),
                                       IAInterval(0.953089922969331, 0.953089922969333)},
                                      {IAInterval(0.239314335249683, 0.239314335249684),
                                       IAInterval(0.230765344947158, 0.230765344947159),
                                       IAInterval(0.769234655052841, 0.769234655052842)}});
  case 10:
  case 11:
    return QintSymT_data<IAInterval>({}, {{IAInterval(0.085662246189585, 0.085662246189586),
                                           IAInterval(0.033765242898423, 0.033765242898424),
                                           IAInterval(0.966234757101575, 0.966234757101577)},
                                          {IAInterval(0.180380786524069, 0.180380786524070),
                                           IAInterval(0.169395306766867, 0.169395306766868),
                                           IAInterval(0.830604693233132, 0.830604693233133)},
                                          {IAInterval(0.233956967286345, 0.233956967286346),
                                           IAInterval(0.380690406958401, 0.380690406958402),
                                           IAInterval(0.619309593041598, 0.619309593041599)}});
  case 12:
  case 13:
    return QintSymT_data<IAInterval>({{IAInterval(0.208979591836734, 0.208979591836735)}},
                                     {{IAInterval(0.064742483084434, 0.064742483084435),
                                       IAInterval(0.025446043828620, 0.025446043828621),
                                       IAInterval(0.974553956171379, 0.974553956171380)},
                                      {IAInterval(0.139852695744638, 0.139852695744639),
                                       IAInterval(0.129234407200302, 0.129234407200303),
                                       IAInterval(0.870765592799697, 0.870765592799698)},
                                      {IAInterval(0.190915025252559, 0.190915025252560),
                                       IAInterval(0.297077424311301, 0.297077424311302),
                                       IAInterval(0.702922575688698, 0.702922575688699)}});
  case 14:
  case 15:
    return QintSymT_data<IAInterval>({}, {{IAInterval(0.050614268145188, 0.050614268145189),
                                           IAInterval(0.019855071751231, 0.019855071751232),
                                           IAInterval(0.980144928248768, 0.980144928248769)},
                                          {IAInterval(0.111190517226687, 0.111190517226688),
                                           IAInterval(0.101666761293186, 0.101666761293187),
                                           IAInterval(0.898333238706813, 0.898333238706814)},
                                          {IAInterval(0.156853322938943, 0.156853322938944),
                                           IAInterval(0.237233795041835, 0.237233795041836),
                                           IAInterval(0.762766204958164, 0.762766204958165)},
                                          {IAInterval(0.181341891689180, 0.181341891689182),
                                           IAInterval(0.408282678752175, 0.408282678752176),
                                           IAInterval(0.591717321247824, 0.591717321247825)}});
  case 16:
  case 17:
    return QintSymT_data<IAInterval>({{IAInterval(0.165119677500629, 0.165119677500630)}},
                                     {{IAInterval(0.040637194180787, 0.040637194180788),
                                       IAInterval(0.015919880246186, 0.015919880246187),
                                       IAInterval(0.984080119753812, 0.984080119753814)},
                                      {IAInterval(0.090324080347428, 0.090324080347429),
                                       IAInterval(0.081984446336682, 0.081984446336683),
                                       IAInterval(0.918015553663317, 0.918015553663318)},
                                      {IAInterval(0.130305348201467, 0.130305348201468),
                                       IAInterval(0.193314283649704, 0.193314283649705),
                                       IAInterval(0.806685716350295, 0.806685716350296)},
                                      {IAInterval(0.156173538520001, 0.156173538520002),
                                       IAInterval(0.337873288298095, 0.337873288298096),
                                       IAInterval(0.662126711701904, 0.662126711701905)}});
  case 18:
  case 19:
    return QintSymT_data<IAInterval>({}, {{IAInterval(0.033335672154344, 0.033335672154345),
                                           IAInterval(0.013046735741414, 0.013046735741415),
                                           IAInterval(0.986953264258585, 0.986953264258586)},
                                          {IAInterval(0.074725674575290, 0.074725674575291),
                                           IAInterval(0.067468316655507, 0.067468316655508),
                                           IAInterval(0.932531683344492, 0.932531683344493)},
                                          {IAInterval(0.109543181257991, 0.109543181257992),
                                           IAInterval(0.160295215850487, 0.160295215850488),
                                           IAInterval(0.839704784149512, 0.839704784149513)},
                                          {IAInterval(0.134633359654998, 0.134633359654999),
                                           IAInterval(0.283302302935376, 0.283302302935377),
                                           IAInterval(0.716697697064623, 0.716697697064624)},
                                          {IAInterval(0.147762112357376, 0.147762112357377),
                                           IAInterval(0.425562830509184, 0.425562830509185),
                                           IAInterval(0.574437169490815, 0.574437169490816)}});
  default:
    Exit("IAQuadrature exact up to degree " + std::to_string(exactUpTo)
         + " not implemented for intervals!");
  }
}

#endif

#endif // QUADRATUREDATAINTERVAL_HPP
