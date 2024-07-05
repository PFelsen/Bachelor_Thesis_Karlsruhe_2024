#include "ParameterImage.hpp"
#include "TimeDate.hpp"
#include <iostream>
#include <cfenv>

bool readPictureLode(const std::string &filename,
                     std::vector<PixelValue> &imageData,
                     unsigned int &width,
                     unsigned int &height) {
  std::vector<ChannelValue> image;
  bool grey =true ;
  if(filename.find("RGBA") != std::string::npos){
    grey = false;
  }
  unsigned error = lodepng::decode(image, width, height, filename);
  if (error) {
    std::cout << "decoder error " << error << ": " << lodepng_error_text(error)
              << "fn: " << filename << std::endl;
    THROW("Could not load Image");
  }
  imageData.clear();
  imageData.resize(width * height);
  string  s = grey ? " Reading Grey" : " Reading RGBA";
  mout << "width: " << width << " height " << height << s << endl;
  mout.StartBlock("ConvertImage");
  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      if (grey) {      // Data is stored as RGBA, for greyscale R=G=B => only take R
        imageData[i * height + j] = image[4 * (i + j * width)];
      }else{// Add all channels
        imageData[i * height + j] = image[4*(i + j * width)];
        imageData[i * height + j] += image[4*(i + j * width) + 1];
        imageData[i * height + j] += image[4*(i + j * width) + 2];
        imageData[i * height + j] += image[4*(i + j * width) + 3];
      }
    }
  }
  mout.EndBlock();
  return grey;
}

void ParameterImage::readAndBroadcastAndCollect(const string &fn_image) {
  //Date start;
  ExchangeBuffer buffer;
  if (PPM->master()) {
    bool grey = readPictureLode(fn_image, imageData, width, height);
    if (grey){
      Dof = 255.0;
    }else{
      Dof = 1020.0;
    }
    PPM->BroadcastDouble(Dof);
    PPM->BroadcastInt(width);
    PPM->BroadcastInt(height);
    PPM->Broadcast(imageData.data(), (width * height) * sizeof(PixelValue), 0);
  }else {
    Dof = PPM->BroadcastDouble();
    width = PPM->BroadcastInt();
    height = PPM->BroadcastInt();
    imageData.resize(width * height);
    PPM->Broadcast(imageData.data(), imageData.size() * sizeof(PixelValue), 0);
  }

  // Time t = Date() - start;
  // pout << "Communicate image: " << t << endl;
}

double ParameterImage::ValueAtPoint(const Point &P) const {
  const auto [low, high] = MinMax;
  std::fesetround( FE_DOWNWARD );
  int ix = toImgCoord(P[0], LeftRight, width);
  int iy = toImgCoord(P[1], BottomTop, height);
  std::fesetround( FE_UPWARD );

  //mout << DOUT(ix) << DOUT(iy) << endl;
  return low + imageData[ix * height + iy] * (high - low) / Dof;
}

int ParameterImage::toImgCoord(double c,
                               const std::pair<double, double> &domainLeftRight,
                               unsigned int pixelNumber) const {
  const auto [domainLeft, domainRight] = domainLeftRight;
  const double size = domainRight - domainLeft;
  if (problemLevel < 0) {
    return int(round((c - domainLeft) / size * pixelNumber));
  } else {
    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    double pow2L = pow(2, problemLevel);
    double pow2L1 = pow(0.5, problemLevel + 1);
    for (int k = 0; k <= pow2L; k++) {
      double fac = (1 + 2 * k) * pow2L1;
      if (abs(dezimal - fac) < pow2L1) {
        //mout << DOUT(ganzzahl) << DOUT(dezimal) << DOUT(pow2L) << DOUT(pow2L1)
        //     << DOUT(k) << DOUT(fac) << DOUT(size) << DOUT(domainRight) << endl;
        return int(round(((ganzzahl + fac) - domainLeft) / size * pixelNumber));
      }
    }
  }
  THROW("corresponding pixel not found")
}
