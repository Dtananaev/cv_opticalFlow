/*
* gaussf.cpp
*
*  Created on: November 23, 2014
*      Author: Denis Tananaev
*/

#include "CMatrix.h"
#include <stdlib.h> 
#include <ctime>
#include <cmath>
#include <string>
#include "CTensor.h"




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

///For colored images
void image2rgb(CTensor<float> image_color, CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue){
	for(int x=0;x<image_color.xSize();x++)
		for(int y=0;y<image_color.ySize();y++)
		{
			red(x,y)=image_color(x,y,0);
			green(x,y)=image_color(x,y,1);
			blue(x,y)=image_color(x,y,2);
		}

}

void rgb2image(CTensor<float>& result_color, CMatrix<float> red, CMatrix<float> green, CMatrix<float> blue){
	for(int x=0;x<result_color.xSize();x++)
		for(int y=0;y<result_color.ySize();y++)
		{
			result_color(x,y,0)=red(x,y);
			result_color(x,y,1)=green(x,y);
			result_color(x,y,2)=blue(x,y);
		}
}

CMatrix<float> Denoise(CMatrix<float> image,int border_size, float alpha, float treshold){
    CMatrix<float> u_k(image);//init u matrix
    CMatrix<float> u_k_new(image);//init u matrix

    //Here implemented denoising with Jacoby method
    /*
      x(k+1)= D^-1 * ( b - M*x(k) )
      b = I -initial image
      D= 1+4*alpha
      M*x(k) = alpha * ( u_k(x-1,y) + u_k(x,y-1) +  u_k(x+1,y) +  u_k(x,y + 1) )

  */

  
    float D = 1+4*alpha;
    float diff; //difference between u(k+1) and u (k) should be less than treshold
    int number_of_pixels = (image.xSize()-2*border_size)*(image.ySize()-2*border_size);; //number of pixels

    do{
         diff = 0; //reset difference

    for(int x=border_size; x<image.xSize()-border_size; x++)
       for(int y=border_size; y<image.ySize()-border_size; y++){
        
        //Energy minimization with Jacoby linear solver
        u_k_new(x,y)=(image(x,y) + alpha * (   u_k(x-1,y) + u_k(x,y-1) +  u_k(x+1,y) +  u_k(x,y + 1) )  )/ D;            
         diff += (u_k_new(x,y) -  u_k(x,y) ) *( u_k_new(x,y) -  u_k(x,y) ); 
    }
     u_k=u_k_new;
     diff= diff/number_of_pixels;
    std::cout<<"diff "<<diff<<"\n";
 }while(diff>treshold);    

 return u_k_new;
}

void applyDenoise( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue,int border_size, float alpha, float treshold){





	red=Neumann_bound(red,border_size);
	green=Neumann_bound(green,border_size);
	blue=Neumann_bound(blue,border_size);
	//Apply denoise


	red=Denoise( red, border_size, alpha, treshold);
	green=Denoise(green, border_size, alpha, treshold);
	blue=Denoise( blue, border_size, alpha, treshold);

	//cut the Neumann boundaries
	red=cut(red, border_size);
	green=cut(green, border_size);
	blue=cut(blue, border_size);
}


int main(int argc, char** argv) {

	std::string fileNameInput;
    int border_size=10;
    float alpha =10;
    float treshold = 0.01;

	if (argc==2){
		fileNameInput=argv[1];

	} else if (argc==3){
		fileNameInput=argv[1];
		alpha=atoi(argv[2]);
	
	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		std::cout<<"Usage: denoiz inputfile <alpha> "<<"\n";
		std::cout<<"The command should contain at least input file name. The default alpha = 10 (smoothing parameter)."<<"\n";
		return 0;    
	}

	std::cout<<"file name is "<<fileNameInput.c_str()<<"\n";

	//Color filter

	CTensor<float> image_color;

	image_color.readFromPPM((fileNameInput+".ppm").c_str());

	CMatrix<float> red(image_color.xSize(),image_color.ySize()),filter_red(image_color.xSize(),image_color.ySize());
	CMatrix<float> green(image_color.xSize(),image_color.ySize()),filter_green(image_color.xSize(),image_color.ySize());
	CMatrix<float> blue(image_color.xSize(),image_color.ySize()),filter_blue(image_color.xSize(),image_color.ySize());

	CTensor<float> result_color(image_color.xSize(),image_color.ySize(),image_color.zSize(),0);

	image2rgb(image_color,red, green, blue);
    
  

	applyDenoise(red, green, blue, border_size, alpha, treshold);

	rgb2image(result_color,red, green, blue);

	result_color.writeToPPM((fileNameInput+"_denoiz.ppm").c_str());



	return 0;
}
