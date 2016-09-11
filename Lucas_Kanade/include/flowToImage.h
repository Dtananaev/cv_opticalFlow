/**
 * @author Nikolaus Mayer, 2015 (mayern@cs.uni-freiburg.de)
 *
 * @brief Image Processing and Computer Graphics
 *        Winter Term 2015/2016
 *        Exercise Sheet 2 (Image Processing part)
 */

#ifndef FLOWTOIMAGE_H__
#define FLOWTOIMAGE_H__

#include "CMatrix.h"
#include "CTensor.h"

/**
 * @brief Convert a 2d vector into an RGB color representation. The
 *        hue encodes the vector's angle while the brightness encodes
 *        the vector's magnitude.
 *
 * @param x Horizontal component of input flow
 * @param y Vertical component of input flow
 * @param R (output parameter) Red output channel, range [0,255]
 * @param G (output parameter) Green output channel, range [0,255]
 * @param B (output parameter) Blue output channel, range [0,255]
 */
void cartesianToRGB(float x, float y, float& R, float& G, float& B);


/**
 * @brief Generate an RGB image as visual representation of a 
 *        two-channel 2d "flow image"
 *
 * @param aFlow Input flow, a tensor with two layers
 * @param aImage (output parameter) RGB flow visualization
 */
void flowToImage(const CTensor<float>& inFlow, CTensor<float>& aImage);


#endif  // FLOWTOIMAGE_H__

