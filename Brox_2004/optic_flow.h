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


class OF{

public:
    OF();
    virtual ~OF();

    //Boundary conditions
    CMatrix<float> Dirichlet_bound(CMatrix<float> aImage,int border_size);
    CMatrix<float> Neumann_bound(CMatrix<float> aImage,int border_size);
    CMatrix<float> cut(CMatrix<float>& image,int border_size);

    //Gaussian Kernel and Filter
    CMatrix<float> Gauss(int sigma);//Kernel
    CMatrix<float> Gfilter(CMatrix<float> Gauss, CMatrix<float> boundary_Image,int border); //Filter

    //Derivatives
    void diffXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy);//  diff -0,5 0 0,5
    void difForwXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy);//  diff -1 1

    //Optic flow estiomation
    void calculateTV(CMatrix<float> u_k, CMatrix<float> v_k, CMatrix<float>& g_1,CMatrix<float>& g_2,CMatrix<float>& g_3,CMatrix<float>& g_4, int border );

    CTensor<float> GaussSeidelHSTV(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold, float gamma);
    CTensor<float> Horn_SchunkOptFlow(CMatrix<float> image1, CMatrix<float> image2, int sigma, bool presmoothing, float alpha, float treshold, float gamma);

private:

};





#endif /* OPTIC_FLOW_H_ */

