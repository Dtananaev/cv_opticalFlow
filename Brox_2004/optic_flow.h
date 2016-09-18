/*
 * File: optic_flow.h
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */



#ifndef OPTIC_FLOW_H_
#define OPTIC_FLOW_H_

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include "CMatrix.h"
#include "CFilter.h" // Csmooth and Cderivative
#include "CTensor.h"
#include "load_sequence.h"
#include "flowToImage.h"
#include "string"
#include "boundary.h"
#include "Gauss_filter.h"

class OF: boundary  {

public:
    OF();
    virtual ~OF();
    //Derivatives
    void diffXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy);//  diff -0,5 0 0,5
    void difForwXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy);//  diff -1 1

    //Optic flow estimation
    void calculateTV(CMatrix<float> u_k, CMatrix<float> v_k, CMatrix<float>& g_1,CMatrix<float>& g_2,CMatrix<float>& g_3,CMatrix<float>& g_4, int border );

    CTensor<float> GaussSeidelHSTV(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold, float gamma); //solver

    CTensor<float> Horn_SchunkOptFlow(CMatrix<float> image1, CMatrix<float> image2, int sigma, bool presmoothing, float alpha, float treshold, float gamma);// whole algorithm

private:
    GaussF filter;
};





#endif /* OPTIC_FLOW_H_ */

