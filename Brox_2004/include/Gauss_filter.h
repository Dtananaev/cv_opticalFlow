/*
 * File: GaussF.cpp
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */

#ifndef GAUSS_FILTER_H_
#define GAUSS_FILTER_H_

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include "CMatrix.h"
#include "CFilter.h" // Csmooth and Cderivative
#include "CTensor.h"
#include "boundary.h"
class GaussF: boundary {

public:
    GaussF();
    virtual ~GaussF();
    
    //Gaussian Kernel and Filter
    CMatrix<float> Gauss(int sigma);//Kernel
    CMatrix<float> Gfilter(CMatrix<float> Gauss, CMatrix<float> boundary_Image,int border); //Filter
private:
};

#endif /* GAUSS_FILTER_H_ */
