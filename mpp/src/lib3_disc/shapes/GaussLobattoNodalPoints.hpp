#ifndef GAUSSLOBATTONODALPOINTS_HPP
#define GAUSSLOBATTONODALPOINTS_HPP

#include "Cell.hpp"
#include "Quadrature.hpp"

const std::unordered_map<int, vector<double>> GAUSSLOBATTO_NODALPOINTS =
    {{1, {0.5}},
     {2, {0.0, 1.0}},
     {3, {0.0, 0.5, 1.0}},
     {4, {0.0, 0.27639320225002103036, 0.72360679774997896964, 1.0}},
     {5, {0.0, 0.17267316464601142810, 0.5, 0.82732683535398857190, 1.0}},
     {6,
      {0.0, 0.11747233803526765357, 0.35738424175967745184, 0.64261575824032254816,
       0.88252766196473234643, 1.0}},
     {7,
      {0.0, 0.084888051860716535064, 0.26557560326464289310, 0.5, 0.73442439673535710690,
       0.91511194813928346494, 1.0}},
     {8,
      {0.0, 0.064129925745196692331, 0.20414990928342884893, 0.39535039104876056562,
       0.60464960895123943438, 0.79585009071657115107, 0.93587007425480330767, 1.0}},
     {9,
      {0.0, 0.050121002294269921344, 0.16140686024463112328, 0.31844126808691092064, 0.5,
       0.68155873191308907936, 0.83859313975536887672, 0.94987899770573007866, 1.0}},
     {10,
      {0.0, 0.040233045916770593086, 0.13061306744724746250, 0.26103752509477775217,
       0.41736052116680648769, 0.58263947883319351231, 0.73896247490522224783,
       0.86938693255275253750, 0.95976695408322940691, 1.0}},
     {11,
      {0.0, 0.032999284795970432834, 0.10775826316842779069, 0.21738233650189749676,
       0.35212093220653030428, 0.5, 0.64787906779346969572, 0.78261766349810250324,
       0.89224173683157220931, 0.96700071520402956717, 1.0}},
     {12,
      {0.0, 0.027550363888558888296, 0.090360339177996660826, 0.18356192348406966117,
       0.30023452951732553387, 0.43172353357253622257, 0.56827646642746377743,
       0.69976547048267446613, 0.81643807651593033883, 0.90963966082200333917,
       0.97244963611144111170, 1.0}},
     {13,
      {0.0, 0.023345076678918044052, 0.076826217674063841567, 0.15690576545912128696,
       0.25854508945433189913, 0.37535653494688000372, 0.5, 0.62464346505311999628,
       0.74145491054566810087, 0.84309423454087871304, 0.92317378232593615843,
       0.97665492332108195595, 1.0}},
     {14,
      {0.0, 0.020032477366369549322, 0.066099473084826374500, 0.13556570045433692971,
       0.22468029853567647234, 0.32863799332864357748, 0.44183406555814806617,
       0.55816593444185193383, 0.67136200667135642252, 0.77531970146432352766,
       0.86443429954566307029, 0.93390052691517362550, 0.97996752263363045068, 1.0}},
     {15,
      {0.0, 0.017377036748080713602, 0.057458977888511850587, 0.11824015502409239965,
       0.19687339726507714444, 0.28968097264316375954, 0.39232302231810288089, 0.5,
       0.60767697768189711911, 0.71031902735683624046, 0.80312660273492285556,
       0.88175984497590760035, 0.94254102211148814941, 0.98262296325191928640, 1.0}},
     {16,
      {0.0, 0.015215976864891033524, 0.050399733453263953503, 0.10399585406909246803,
       0.17380564855875345527, 0.25697028905643119411, 0.35008476554961839595,
       0.44933686323902527608, 0.55066313676097472392, 0.64991523445038160405,
       0.74302971094356880589, 0.82619435144124654473, 0.89600414593090753197,
       0.94960026654673604650, 0.98478402313510896648, 1.0}},
     {17,
      {0.0, 0.013433911684290842922, 0.044560002042213202188, 0.092151874389114846447,
       0.15448550968615764730, 0.22930730033494923044, 0.31391278321726147905,
       0.40524401324084130585, 0.5, 0.59475598675915869415, 0.68608721678273852095,
       0.77069269966505076956, 0.84551449031384235270, 0.90784812561088515355,
       0.95543999795778679781, 0.98656608831570915708, 1.0}},
     {18,
      {0.0, 0.011947221293900728568, 0.039675407326233063081, 0.082203232390954893143,
       0.13816033535837865935, 0.20574758284066911941, 0.28279248154393801233,
       0.36681867356085950792, 0.45512545325767394449, 0.54487454674232605551,
       0.63318132643914049208, 0.71720751845606198767, 0.79425241715933088059,
       0.86183966464162134065, 0.91779676760904510686, 0.96032459267376693692,
       0.98805277870609927143, 1.0}},
     {19,
      {0.0, 0.010694116888959952424, 0.035549235923706878141, 0.073769711101676953457,
       0.12425289872369349292, 0.18554593136738975112, 0.25588535715964324861,
       0.33324757608775069485, 0.41540698829535921431, 0.5, 0.58459301170464078569,
       0.66675242391224930515, 0.74411464284035675139, 0.81445406863261024888,
       0.87574710127630650708, 0.92623028889832304654, 0.96445076407629312186,
       0.98930588311104004758, 1.0}}};

template<typename T, int sDim, int tDim>
vector<PointT<T, sDim, tDim>> GaussLobattoNodalPointsInterval(const Cell &c, int npcount) {
  if (npcount == 1) { return vector<PointT<T, sDim, tDim>>{(c[0] + c[1]) / 2.0}; }
  vector<double> Q = GAUSSLOBATTO_NODALPOINTS.at(npcount);
  vector<PointT<T, sDim, tDim>> result;
  for (const Point &np : Q) {
    PointT<T, sDim, tDim> p = c[0] + (c[1] - c[0]) * np[0];
    result.push_back(p);
  }
  return result;
}

template<typename T, int sDim, int tDim>
vector<PointT<T, sDim, tDim>> GaussLobattoNodalPointsQuadrilateral(const Cell &c, int npcount) {
  vector<double> Q = GAUSSLOBATTO_NODALPOINTS.at(npcount);
  vector<PointT<T, sDim, tDim>> result;
  for (const double &np_i : Q) {
    for (const double &np_j : Q) {
      PointT<T, sDim, tDim> v = c[0] + np_j * (c[1] - c[0]);
      PointT<T, sDim, tDim> w = c[3] + np_j * (c[2] - c[3]);
      result.push_back(v + np_i * (w - v));
    }
  }
  return result;
}

template<typename T, int sDim, int tDim>
vector<PointT<T, sDim, tDim>> GaussLobattoNodalPointsHexaedron(const Cell &c, int npcount) {
  if (npcount < 2) { THROW("GaussLobatto need at leat 2 Nodalpoints.") }
  vector<double> Q = GAUSSLOBATTO_NODALPOINTS.at(npcount);
  vector<PointT<T, sDim, tDim>> result;
  for (const double &np_i : Q) {
    for (const double &np_j : Q) {
      for (const double &np_k : Q) {
        PointT<T, sDim, tDim> v = c[0] + np_k * (c[1] - c[0]);
        PointT<T, sDim, tDim> w = c[3] + np_k * (c[2] - c[3]);
        PointT<T, sDim, tDim> x = c[4] + np_k * (c[5] - c[4]);
        PointT<T, sDim, tDim> y = c[7] + np_k * (c[6] - c[7]);
        PointT<T, sDim, tDim> r = v + np_j * (w - v);
        PointT<T, sDim, tDim> s = x + np_j * (y - x);
        result.push_back(r + np_i * (s - r));
      }
    }
  }
  return result;
}

#endif // GAUSSLOBATTONODALPOINTS_HPP