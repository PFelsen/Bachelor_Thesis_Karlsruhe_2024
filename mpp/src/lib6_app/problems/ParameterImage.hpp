#ifndef PARAMETERIMAGE_HPP
#define PARAMETERIMAGE_HPP

#include "lodepng.hpp"
#include "GlobalDefinitions.hpp"
#include "ExchangeBuffer.hpp"
#include "Parallel.hpp"
#include "Config.hpp"
#include "Point.hpp"

using ChannelValue = unsigned char;
using PixelValue = unsigned short;

bool readPictureLode(const std::string &filename,
                     std::vector<PixelValue> &imageData,
                     unsigned &width,
                     unsigned &height);

class ParameterImage {
private:

  vector<PixelValue> imageData;
  unsigned width = 0;
  unsigned height = 0;
  int problemLevel = -1;
  double Dof = -1;
  void readAndBroadcastAndCollect(const std::string &fn_image);

  int toImgCoord(double c, const std::pair<double, double> &CoordPair, unsigned int pixelNumber) const;

  std::pair<double, double> MinMax;
  std::pair<double, double> LeftRight;
  std::pair<double, double> BottomTop;
public:
  ParameterImage(const std::string &parameterName,
                 const std::string &fn_image,
                 std::pair<double, double> LeftRight,
                 std::pair<double, double> BottomTop)
      : LeftRight(LeftRight), BottomTop(BottomTop) {
    readAndBroadcastAndCollect(fn_image);
    if (LeftRight.first > LeftRight.second || (BottomTop.first > BottomTop.second)) {
      THROW("error")
    }
    Config::Get("Min" + parameterName, MinMax.first);
    Config::Get("Max" + parameterName, MinMax.second);
    Config::Get("ProblemLevel", problemLevel);
    mout << parameterName + " model read   "
         << fn_image << "' (size: " << width << "x"
         << height << ")" << endl;
  }
  double getDof(){
    return Dof;
  }
  double ValueAtPoint(const Point &P) const;
};

#endif

