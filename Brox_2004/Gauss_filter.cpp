/*
 * File: GaussF.cpp
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */


#include "Gauss_filter.h"


GaussF::GaussF(){
}
GaussF::~GaussF(){
}


//Gaussian Kernel 2D witch acc 3 sigma
CMatrix<float> GaussF::Gauss(int sigma){

	size_t filter_size = size_t(6*sigma+1);
	std::cout<<"filter_size_Gauss= "<<filter_size<<" pixels"<<std::endl;
	CMatrix<float> filter(filter_size,filter_size);
	if ( filter_size % 2 == 0 )
		++filter_size;
	int m=(filter_size-1)/2;
	filter.fill(0);
	float n=0;
	float sum=0;
	float sum1=0;
	n=2*sigma*sigma;
	for(int x=-m;x<=m;x++)
		for(int y=-m;y<=m;y++)
		{
			filter(x+m,y+m)=(exp(-((x*x+y*y)/n)))/n*3.1415926536;
			sum+=filter(x+m,y+m);

		}	

		// normalize the Kernel
		for(size_t i = 0; i < filter_size; ++i)
			for(size_t j = 0; j < filter_size; ++j)
			{  filter(i,j) /= sum;
		sum1+=filter(i,j);

		}
		return filter;
}

//Gaussian filter
CMatrix<float> GaussF::Gfilter(CMatrix<float> Gauss, CMatrix<float> boundary_Image,int border){
	CMatrix<float> image( boundary_Image.xSize(),boundary_Image.ySize());
	image.fill(0);
	//int border=sigma*3;

	//Filtering
	for(int x=border;x<boundary_Image.xSize()-border;x++)
		for(int y=border;y<boundary_Image.ySize()-border;y++)
		{
			for(int i=0;i<Gauss.xSize();i++)
				for (int j=0;j<Gauss.ySize();j++)
				{
					image(x,y)+=boundary_Image(x-border+i,y-border+j)*Gauss(i,j);
				}
		}

		return image;
}
