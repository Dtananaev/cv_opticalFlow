/*
 * File: boundary.h
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */


#ifndef BOUNDARY_H_
#define BOUNDARY_H_ 

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include "CMatrix.h"
#include "CFilter.h" // Csmooth and Cderivative
#include "CTensor.h"



class boundary{
public:
    boundary();
    virtual ~boundary();


    //Boundary conditions
    CMatrix<float> Dirichlet_bound(CMatrix<float> aImage,int border_size);
    CMatrix<float> Neumann_bound(CMatrix<float> aImage,int border_size);
    CMatrix<float> cut(CMatrix<float>& image,int border_size);

private:

};

#endif /* BOUNDARY_H_ */



