/**
 * Image Processing and Computer Graphics
 *
 * Created on:  September 10, 2016
 *
 * File: Lucas_Kanade.cpp
 *
 * Author: Denis Tananaev
 *  
 */

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include "CMatrix.h"
#include "CFilter.h" // Csmooth and Cderivative
#include "CTensor.h"
#include "load_sequence.h"
#include "flowToImage.h"
#include "string"
//Cut Neumann boundaries
CMatrix<float> cut(CMatrix<float>& image,int border_size){ 

	CMatrix<float> realimage(image.xSize()-2*border_size,image.ySize()-2*border_size);
	for(int x=0;x<realimage.xSize();x++)
		for(int y=0; y<realimage.ySize();y++)
		{
			realimage(x,y)=image(x+border_size,y+border_size);
		}

		return realimage;
}

// Neumann boundry conditions
CMatrix<float> Neumann_bound(CMatrix<float> aImage,int border_size){

	CMatrix<float> result(aImage.xSize()+border_size*2,aImage.ySize()+border_size*2);
	result.fill(0);
	//center matrix
	for(int x=0;x<aImage.xSize();x++)
		for(int y=0;y<aImage.ySize();y++)
		{
			result(x+border_size,y+border_size)=aImage(x,y);
		}
		//Top
		for(int x=0;x<aImage.xSize();x++)
			for(int y=0;y<border_size;y++)
			{
				result(x+border_size,y)=aImage(x,border_size-1-y);
			}
			//Bottom
			for(int x=0;x<aImage.xSize();x++)
				for(int y=0;y<border_size;y++)
				{
					result(x+border_size,y+aImage.ySize()+border_size)=aImage(x,aImage.ySize()-1-y);
				}
				//left side
				for(int x=0;x<border_size;x++)
					for(int y=0;y<aImage.ySize();y++)
					{
						result(x,y+border_size)=aImage(border_size-1-x,y);
					}

					//right side
					for(int x=0;x<border_size;x++)
						for(int y=0;y<aImage.ySize();y++)
						{
							result(x+aImage.xSize()+border_size,y+border_size)=aImage(aImage.xSize()-1-x,y);
						}
						//up left square
						for(int x=0;x<border_size;x++)
							for(int y=0;y<border_size;y++)
							{
								result(x,y)=aImage(0,0);
							}
							//up right square
							for(int x=aImage.xSize()-1;x<(aImage.xSize()+border_size);x++)
								for(int y=0;y<border_size;y++)
								{
									result(x+border_size,y)=aImage(aImage.xSize()-1,0);
								}
								//down left square
								for(int x=0;x<border_size;x++)
									for(int y=aImage.ySize()-1;y<(aImage.ySize()+border_size);y++)
									{
										result(x,y+border_size)=aImage(0,aImage.ySize()-1);
									}
									//down right square
									for(int x=aImage.xSize()-1;x<(aImage.xSize()+border_size);x++)
										for(int y=aImage.ySize()-1;y<(aImage.ySize()+border_size);y++)
										{
											result(x+border_size,y+border_size)=aImage(aImage.xSize()-1,aImage.ySize()-1);
										}		
										return result;
}


//Gaussian Kernel 2D witch acc 3 sigma
CMatrix<float> Gauss(int sigma){

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
CMatrix<float> Gfilter(CMatrix<float> Gauss, CMatrix<float> boundary_Image,int border){
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

//  diff -0,5 0 0,5
void diffXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy){
       int w =Image.xSize();
       int h =Image.ySize();
       Image=Neumann_bound(Image,1);
    dx.fill(0);
    dy.fill(0);
    for(int x=0; x<w; x++)
    		for(int y=0; y<h; y++) {
            int a=x+1;
            int b=y+1;
            dx(x,y)=-0.5*Image(a-1,b)+0.5*Image(a+1,b);
            dy(x,y)=-0.5*Image(a,b-1)+0.5*Image(a,b+1);     
 
    }

}

CTensor<float> Lucas_KanadeOptFlow(CMatrix<float> image1, CMatrix<float> image2, int sigma, bool presmoothing){
    
    int width =image1.xSize();
    int height=image1.ySize();

    CTensor<float> result(width,height,2); //u and v optic flow result
    int border_size=3*sigma;
    CMatrix<float> kernel;
 
    if(sigma>0){
      kernel=Gauss(sigma);//create Gauss kernel
        //Pre-smooth  image with Gauss kernel if presmoothing is true
        if(presmoothing==true){
        //apply Neumann boundary condition 
        image1=Neumann_bound(image1,border_size);
	    image2=Neumann_bound(image2,border_size); 
        //Pre-smoothing initial images
        image1=Gfilter(kernel, image1, border_size);
        image2=Gfilter(kernel, image2, border_size);
        image1=cut(image1, border_size);
        image2=cut(image2, border_size);    
        }
    }


    /// Derivative filter [-0.5, 0, 0.5]
    CDerivative<float> derivative(3);

    // apply derivatives
    CMatrix<float> dx(width,height);
   CMatrix<float> dy(width,height); 
   // CMatrix<float> dx(image1);
   // CMatrix<float> dy(image1); 
    //NFilter::filter(dx, derivative, 1);
    //NFilter::filter(dy, 1, derivative);   
    CMatrix<float> dt(image2);  

    diffXY(image1, dx, dy);
std::cout<<"511"<<"\n"; 
   // diffT(image1, image2, dt);
std::cout<<"5"<<"\n"; 
    //Use the Gauss kernel in order to smooth derivatives

    /*
        The Lucas-Kanade equations:

        | Ixx Ixy | |u| = - Ixt
        | Ixy Iyy | |v|   - Iyt


    */
    CMatrix<float> Ixx(width,height);
    CMatrix<float> Ixy(width,height);  
    CMatrix<float> Ixt(width,height);  
    CMatrix<float> Iyy(width,height);
    CMatrix<float> Iyt(width,height);
   for(int x=0; x<width; x++)
    		for(int y=0;y<height;y++) { 
        dt(x,y)-=image1(x,y);
        Ixx(x,y) = dx(x,y)*dx(x,y);
        Ixy(x,y) = dx(x,y)*dy(x,y);
        Iyy(x,y) = dy(x,y)*dy(x,y);
        Ixt(x,y) = dx(x,y)*dt(x,y);
        Iyt(x,y) = dy(x,y)*dt(x,y);
    }

   if(sigma>0){
    
    //apply Neumann boundary condition 
      Ixx=Neumann_bound(Ixx,border_size);
	  Ixy=Neumann_bound(Ixy,border_size); 
      Iyy=Neumann_bound(Iyy,border_size);
	  Ixt=Neumann_bound(Ixt,border_size); 
	  Iyt=Neumann_bound(Iyt,border_size); 
    //apply Gauss smoothing

        Ixx=Gfilter(kernel, Ixx, border_size);
        Ixy=Gfilter(kernel, Ixy, border_size);
        Iyy=Gfilter(kernel, Iyy, border_size);      
        Ixt=Gfilter(kernel, Ixt, border_size);
        Iyt=Gfilter(kernel, Iyt, border_size);

    //cut the Neumann boundaries
	Ixx=cut(Ixx, border_size);
	Ixy=cut(Ixy, border_size);
	Iyy=cut(Iyy, border_size);
	Ixt=cut(Ixt, border_size);
	Iyt=cut(Iyt, border_size);
std::cout<<"7"<<"\n"; 

    }



    //Calculating Lucas-Kanade optic flow
    /*
        u = ( -Ixt + (Iyt*Ixy)/Iyy )/( (Ixy*Ixy)/Iyy + Ixx  )
       v = (-Iyt - Ixy*u ) / Iyy

   */

   float u,v;


    for(int x=0; x<width; x++)
    		for(int y=0;y<height;y++) {  

            u = (-1*Ixt(x,y) + ( Iyt(x,y)*Ixy(x,y)/Iyy(x,y) ) )/( (Ixy(x,y)*Ixy(x,y)/Iyy(x,y)) + Ixx(x,y) );
            v= ( -1* Iyt(x,y) - Ixy(x,y)* u ) / Iyy(x,y);

        //Result
        result(x,y,0)=u;
        result(x,y,1)=v;

    }

    return result; 
}


int main(int argc, char** argv) {  


    std::string folderNameInput;
	int sigma=5;



	if (argc==2){
		folderNameInput=argv[1];

	} else if (argc==3){
		folderNameInput=argv[1];
		sigma=atoi(argv[2]);

	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		std::cout<<"Usage: Lucas_Kanade inputfolder <variance of Gauss kernell>"<<"\n";
		std::cout<<"The command should contain at least input file name. The default Gauss kernel sigma=1."<<"\n";
		return 0;    
	}


    CVector<CMatrix<float> > seq;
    
   // seq = loadSequence("resources/cropped-street/t.txt");
    seq = loadSequence("resources/yos/t.txt");
  // seq = loadSequence("resources/gsalesman/t.txt");
      CMatrix<float> img1;
      CMatrix<float> img2;
      bool presmoothing=false;
       CTensor<float> opticFlow;

    for (int i = 0; i < seq.size()-1; ++i){

	    img1 = seq(i);
	    img2 = seq(i+1);
        opticFlow(img1.xSize(),img1.ySize(),2);

        opticFlow=Lucas_KanadeOptFlow(img1 ,img2, sigma,  presmoothing);

    CTensor<float> LucasKanadeFlowRGB(img1.xSize(), img1.ySize(),3);

    flowToImage(opticFlow, LucasKanadeFlowRGB);

	LucasKanadeFlowRGB.writeToPPM(("result/"+folderNameInput+ std::to_string(i) + ".ppm").c_str());

    }

  return 0;
}

