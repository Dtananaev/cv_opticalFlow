/**
 * @author Nikolaus Mayer, 2015 (mayern@cs.uni-freiburg.de)
 *
 * @brief Image Processing and Computer Graphics
 *        Winter Term 2015/2016
 *        Exercise Sheet 2 (Image Processing part)
 */

/// System/STL
#include <cmath>
/// Local files
#include "CMatrix.h"
#include "CTensor.h"
#include "flowToImage.h"


#define FACTOR 1.f
#define PI 3.1415926536f



// cartesianToRGB
void cartesianToRGB (float x, float y, float& R, float& G, float& B) 
{
  /// Euclidean length of the flow vector
  float radius = std::sqrt(x * x + y * y);
  /// Limit radius to 1.0
  if (radius > 1.f) radius = 1.f;

  /// Compute flow angle
  float phi;
  if (x == 0.0f)
    if (y >= 0.0f) phi = 0.5f * PI;
    else phi = 1.5f * PI;
  else if (x > 0.0f)
    if (y >= 0.0f) phi = std::atan(y/x);
    else phi = 2.0f * PI + std::atan(y/x);
  else phi = PI + std::atan(y/x);

  // Weights for linear interpolation
  float alpha, beta;    
  phi *= 0.5f;
  /// Interpolation between red (0) and blue (0.25 * PI)
  if ((phi >= 0.0f) && (phi < 0.125f * PI)) {
    beta  = phi / (0.125f * PI);
    alpha = 1.0f - beta;
    R = (int)(radius * (alpha * 255.0f + beta * 255.0f));
    G = (int)(radius * (alpha *   0.0f + beta *   0.0f));
    B = (int)(radius * (alpha *   0.0f + beta * 255.0f));
  }
  if ((phi >= 0.125f * PI) && (phi < 0.25f * PI)) {
    beta  = (phi-0.125f * PI) / (0.125f * PI);
    alpha = 1.0f - beta;
    R = (int)(radius * (alpha * 255.0f + beta *  64.0f));
    G = (int)(radius * (alpha *   0.0f + beta *  64.0f));
    B = (int)(radius * (alpha * 255.0f + beta * 255.0f));
  }
  /// Interpolation between blue (0.25 * PI) and green (0.5 * PI)
  if ((phi >= 0.25f * PI) && (phi < 0.375f * PI)) {
    beta  = (phi - 0.25f * PI) / (0.125f * PI);
    alpha = 1.0f - beta;
    R = (int)(radius * (alpha *  64.0f + beta *   0.0f));
    G = (int)(radius * (alpha *  64.0f + beta * 255.0f));
    B = (int)(radius * (alpha * 255.0f + beta * 255.0f));
  }
  if ((phi >= 0.375f * PI) && (phi < 0.5f * PI)) {
    beta  = (phi - 0.375f * PI) / (0.125f * PI);
    alpha = 1.0f - beta;
    R = (int)(radius * (alpha *   0.0f + beta *   0.0f));
    G = (int)(radius * (alpha * 255.0f + beta * 255.0f));
    B = (int)(radius * (alpha * 255.0f + beta *   0.0f));
  }
  /// Interpolation between green (0.5 * PI) and yellow (0.75 * PI)
  if ((phi >= 0.5f * PI) && (phi < 0.75f * PI)) {
    beta  = (phi - 0.5f * PI) / (0.25f * PI);
    alpha = 1.0f - beta;
    R = (int)(radius * (alpha *   0.0f + beta * 255.0f));
    G = (int)(radius * (alpha * 255.0f + beta * 255.0f));
    B = (int)(radius * (alpha *   0.0f + beta *   0.0f));
  }
  /// Interpolation between yellow (0.75 * PI) and red (Pi)
  if ((phi >= 0.75f * PI) && (phi <= PI)) {
    beta  = (phi - 0.75f * PI) / (0.25f * PI);
    alpha = 1.0f - beta;
    R = (int)(radius * (alpha * 255.0f + beta * 255.0f));
    G = (int)(radius * (alpha * 255.0f + beta *   0.0f));
    B = (int)(radius * (alpha *   0.0f + beta *   0.0f));
  }
  if (R < 0.f) R = 0.f;
  if (G < 0.f) G = 0.f;
  if (B < 0.f) B = 0.f;
  if (R > 255.f) R = 255.f;
  if (G > 255.f) G = 255.f;
  if (B > 255.f) B = 255.f;
}

void flowToImage(const CTensor<float>& aFlow, CTensor<float>& aImage) 
{
  /// Scaling factor determined by experiment; prevent over- or
  /// undersaturated image results
  float aFactor = FACTOR*std::sqrt(0.5f)*0.5f;

  int aSize = aFlow.xSize()*aFlow.ySize();
  float R, G, B;
  /// Convert each pixel's flow to RGB color
  for (int i = 0; i < aSize; i++) {
    R = G = B = 0.f;
    cartesianToRGB(aFactor*aFlow.data()[i],
                   aFactor*aFlow.data()[i+aSize],
                   R, G, B);
    aImage.data()[i]         = R;
    aImage.data()[i+aSize]   = G;
    aImage.data()[i+2*aSize] = B;
  }
}
