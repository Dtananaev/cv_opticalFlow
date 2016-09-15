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


enum method{
    Jacoby=1,
    Gauss_Seidel=2,
    SOR=3

};

//Dirichlet boundary conditions
CMatrix<float> Dirichlet_bound(CMatrix<float> aImage,int border_size){

	CMatrix<float> result(aImage.xSize()+border_size*2,aImage.ySize()+border_size*2);
	result.fill(0);
	//center matrix
	for(int x=0;x<aImage.xSize();x++)
		for(int y=0;y<aImage.ySize();y++)
		{
			result(x+border_size,y+border_size)=aImage(x,y);
		}

return result;

}
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

//  diff -1 1
void difForwXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy){
       int w =Image.xSize();
       int h =Image.ySize();
       Image=Neumann_bound(Image,1);
    dx.fill(0);
    dy.fill(0);
    for(int x=0; x<w; x++)
    		for(int y=0; y<h; y++) {
            int a=x+1;
            int b=y+1;
            dx(x,y)=-1*Image(a,b)+1*Image(a+1,b);
            dy(x,y)=-1*Image(a,b)+1*Image(a,b+1);     
 
    }

}


//Jacoby method
CTensor<float> JacobyHS(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold){

         /*Jacoby method
            Linear system:
            Ax=b;  Split A=D+M; D-diagonal part of matrix ; M- outer diagonal part
            Iterative equation:
            x(k+1)=D^-1 * (b-Mx(k)) 
        */
   
    int width =image1.xSize();
    int height=image1.ySize();
    CTensor<float> result(width,height,2); //u and v optic flow result

    int number_of_pixels=width*height;
  // apply derivatives
    CMatrix<float> Ix(width,height);
    CMatrix<float> Iy(width,height); 
    diffXY(image1, Ix, Iy);


    CMatrix<float> Iz(image2);  
      for(int x=0; x<width; x++)
    	for(int y=0;y<height;y++) { 
        Iz(x,y)-=image1(x,y);
}
   
    //To make computation simpler  and have all matrices of the same size apply Neumann boundary conditions with border size 1 to derivatives and image1
    image1=Dirichlet_bound(image1,1);
    Iz=Dirichlet_bound(Iz,1);
    Ix=Dirichlet_bound(Ix,1);
    Iy=Dirichlet_bound(Iy,1);    
    //Horn-Schunck optic flow with Jacoby method
   //  CMatrix<float> u_k(image1);
    // CMatrix<float> v_k(image1); 
     CMatrix<float> u_k(image1.xSize(),image1.ySize(),1);
     CMatrix<float> v_k(image1.xSize(),image1.ySize(),1);        
     CMatrix<float> u_k_new(image1.xSize(),image1.ySize(),0);
     CMatrix<float> v_k_new(image1.xSize(),image1.ySize(),0);
     float diff_u,diff_v;


   
                    
  do{ 
       diff_u=0;
       diff_v=0;  
    for(int x=1; x<image1.xSize()-1; x++)
    		for(int y=1;y<image1.ySize()-1;y++) {  
//Gauss-Seidel method

                //Lower triangle (x-1,y) and (x,y-1)

            u_k_new(x,y)= (1 / (4*alpha + Ix(x,y)*Ix(x,y)) ) * ( alpha* (u_k(x-1,y)+u_k(x+1,y)+u_k(x,y-1)+u_k(x,y+1))  -Ix(x,y)*Iy(x,y)*v_k(x,y) - Ix(x,y)*Iz(x,y) );
            

            v_k_new(x,y) = (1/ (4*alpha+ Iy(x,y)*Iy(x,y)) )*(alpha* (v_k(x-1,y)+v_k(x+1,y)+v_k(x,y-1)+v_k(x,y+1)) - Iy(x,y)*Ix(x,y)*u_k(x,y) -Iy(x,y)*Iz(x,y) );

        diff_u+= (u_k_new(x,y)-u_k(x,y))*(u_k_new(x,y)-u_k(x,y));
        diff_v+= (v_k_new(x,y)-v_k(x,y))*(v_k_new(x,y)-v_k(x,y));    
            
          

    }
         u_k=u_k_new;
         v_k=v_k_new;  
          diff_u=  diff_u/ number_of_pixels;
          diff_v=  diff_v/ number_of_pixels;
         std::cout<<"diff_u "<<diff_u<<"\n";
         std::cout<<"diff_v "<<diff_v<<"\n";
} while (diff_u >treshold && diff_v> treshold);

        u_k_new=cut(u_k_new,1);
        v_k_new=cut(v_k_new,1); 

       result.putMatrix(u_k_new,0);      
       result.putMatrix(v_k_new,1);      

    return result; 
}




//Gauss-Seidel method
CTensor<float> GaussSeidelHS(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold){


         /*Gauss-Seidel method
            Linear system:
            Ax=b;  Split A=D+L+U; D-diagonal part of matrix ; L- lower diagonal part; U - upper diagonal part
            Iterative equation:
            x(k+1)=D^-1 * (b-L*x(k+1)+U*x(k)) 
        */
    
    
    int width =image1.xSize();
    int height=image1.ySize();
    CTensor<float> result(width,height,2); //u and v optic flow result

    int number_of_pixels=width*height;
  // apply derivatives
    CMatrix<float> Ix(width,height);
    CMatrix<float> Iy(width,height); 
    diffXY(image1, Ix, Iy);


    CMatrix<float> Iz(image2);  
      for(int x=0; x<width; x++)
    	for(int y=0;y<height;y++) { 
        Iz(x,y)-=image1(x,y);
}
   
    //To make computation simpler  and have all matrices of the same size apply Neumann boundary conditions with border size 1 to derivatives and image1
    image1=Dirichlet_bound(image1,1);
    Iz=Dirichlet_bound(Iz,1);
    Ix=Dirichlet_bound(Ix,1);
    Iy=Dirichlet_bound(Iy,1);    
    //Horn-Schunck optic flow with Jacoby method
     //CMatrix<float> u_k(image1);
     //CMatrix<float> v_k(image1);
     CMatrix<float> u_k(image1.xSize(),image1.ySize(),1);
     CMatrix<float> v_k(image1.xSize(),image1.ySize(),1);      
     CMatrix<float> u_k_new(image1.xSize(),image1.ySize(),0);
     CMatrix<float> v_k_new(image1.xSize(),image1.ySize(),0);
     float diff_u,diff_v;


   
                    
  do{ 
       diff_u=0;
       diff_v=0;  
    for(int x=1; x<image1.xSize()-1; x++)
    		for(int y=1;y<image1.ySize()-1;y++) {  
//Gauss-Seidel method

                //Lower triangle (x-1,y) and (x,y-1)

            u_k_new(x,y)= (1 / (4*alpha + Ix(x,y)*Ix(x,y)) ) * ( alpha* (u_k_new(x-1,y)+u_k(x+1,y)+u_k_new(x,y-1)+u_k(x,y+1))  -Ix(x,y)*Iy(x,y)*v_k(x,y) - Ix(x,y)*Iz(x,y) );
            

            v_k_new(x,y) = (1/ (4*alpha+ Iy(x,y)*Iy(x,y)) )*(alpha* (v_k_new(x-1,y)+v_k(x+1,y)+v_k_new(x,y-1)+v_k(x,y+1)) - Iy(x,y)*Ix(x,y)*u_k_new(x,y) -Iy(x,y)*Iz(x,y) );

        diff_u+= (u_k_new(x,y)-u_k(x,y))*(u_k_new(x,y)-u_k(x,y));
        diff_v+= (v_k_new(x,y)-v_k(x,y))*(v_k_new(x,y)-v_k(x,y));    
            
          

    }
         u_k=u_k_new;
         v_k=v_k_new;  
          diff_u=  diff_u/ number_of_pixels;
          diff_v=  diff_v/ number_of_pixels;
         std::cout<<"diff_u "<<diff_u<<"\n";
         std::cout<<"diff_v "<<diff_v<<"\n";
} while (diff_u >treshold && diff_v> treshold);

        u_k_new=cut(u_k_new,1);
        v_k_new=cut(v_k_new,1); 

       result.putMatrix(u_k_new,0);      
       result.putMatrix(v_k_new,1);      

    return result; 

}


//Succesive over-relaxation (SOR)
CTensor<float> SORHS(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold){


         /*Succesive over-relaxation method
            Linear system:
            Ax=b;  Split A=D+L+U; D-diagonal part of matrix ; L- lower diagonal part; U - upper diagonal part
            Iterative equation with over-relaxation:
            x(k+1)= (1-w)*x(k)  + w* D^-1 * (b-L*x(k+1)+U*x(k)) 

            w between 0 and 1
        */
    float w = 1.67;
    
    int width =image1.xSize();
    int height=image1.ySize();
    CTensor<float> result(width,height,2); //u and v optic flow result

    int number_of_pixels=width*height;
  // apply derivatives
    CMatrix<float> Ix(width,height);
    CMatrix<float> Iy(width,height); 
    diffXY(image1, Ix, Iy);


    CMatrix<float> Iz(image2);  
      for(int x=0; x<width; x++)
    	for(int y=0;y<height;y++) { 
        Iz(x,y)-=image1(x,y);
}
   
    //To make computation simpler  and have all matrices of the same size apply Neumann boundary conditions with border size 1 to derivatives and image1
    image1=Dirichlet_bound(image1,1);
    Iz=Dirichlet_bound(Iz,1);
    Ix=Dirichlet_bound(Ix,1);
    Iy=Dirichlet_bound(Iy,1);    
    //Horn-Schunck optic flow with Jacoby method
     //CMatrix<float> u_k(image1);
     //CMatrix<float> v_k(image1);  
     CMatrix<float> u_k(image1.xSize(),image1.ySize(),1);
     CMatrix<float> v_k(image1.xSize(),image1.ySize(),1);       
     CMatrix<float> u_k_new(image1.xSize(),image1.ySize(),0);
     CMatrix<float> v_k_new(image1.xSize(),image1.ySize(),0);
     float diff_u,diff_v;


   
                    
  do{ 
       diff_u=0;
       diff_v=0;  
    for(int x=1; x<image1.xSize()-1; x++)
    		for(int y=1;y<image1.ySize()-1;y++) {  
//Gauss-Seidel method

                //Lower triangle (x-1,y) and (x,y-1)

            u_k_new(x,y)=(1-w)*u_k(x,y)+ w * (1 / (4*alpha + Ix(x,y)*Ix(x,y)) ) * ( alpha* (u_k_new(x-1,y)+u_k(x+1,y)+u_k_new(x,y-1)+u_k(x,y+1))  -Ix(x,y)*Iy(x,y)*v_k(x,y) - Ix(x,y)*Iz(x,y) );
            

            v_k_new(x,y) =(1-w)*v_k(x,y)+ w * (1/ (4*alpha+ Iy(x,y)*Iy(x,y)) )*(alpha* (v_k_new(x-1,y)+v_k(x+1,y)+v_k_new(x,y-1)+v_k(x,y+1)) - Iy(x,y)*Ix(x,y)*u_k_new(x,y) -Iy(x,y)*Iz(x,y) );

        diff_u+= (u_k_new(x,y)-u_k(x,y))*(u_k_new(x,y)-u_k(x,y));
        diff_v+= (v_k_new(x,y)-v_k(x,y))*(v_k_new(x,y)-v_k(x,y));    
            
          

    }
         u_k=u_k_new;
         v_k=v_k_new;  
          diff_u=  diff_u/ number_of_pixels;
          diff_v=  diff_v/ number_of_pixels;
            std::cout<<"diff_u "<< diff_u<<"\n";
            std::cout<<"diff_v "<< diff_v<<"\n";
} while (diff_u >treshold && diff_v> treshold);

        u_k_new=cut(u_k_new,1);
        v_k_new=cut(v_k_new,1); 

       result.putMatrix(u_k_new,0);      
       result.putMatrix(v_k_new,1);      

    return result; 

}


CTensor<float> Horn_SchunkOptFlow(CMatrix<float> image1, CMatrix<float> image2, int sigma, bool presmoothing, float alpha, float treshold, int method_choice){
    
    int width =image1.xSize();
    int height=image1.ySize();
   // int number_of_pixels=width*height;

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


    if(method_choice==method::Jacoby){
        result=JacobyHS(image1, image2,  alpha,  treshold);
    }else if(method_choice==method::Gauss_Seidel){
        result=GaussSeidelHS( image1,  image2,  alpha,  treshold);
    }else if(method_choice==method::SOR){
        result = SORHS( image1, image2,  alpha,  treshold);
    }
return result;
  
}


int main(int argc, char** argv) {  


    std::string folderNameInput;

	int sigma=1;



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
    

    int method_choice;

    std::cout<<"Choose solver: [1] - Jacoby; [2] - Gauss-Seidel; [3] - SOR "<<"\n";
    std::cin>>method_choice;
    if(method_choice ==1 || method_choice==2 || method_choice==3 ){
        std::cout<<"Processing..."<<"\n";
    }else{
         std::cout<<"Wrong input! By default activated  Jacoby solver"<<"\n";
         method_choice=method::Jacoby;

    }

    CVector<CMatrix<float> > seq;
    
    //seq = loadSequence("resources/cropped-street/t.txt");
 //   seq = loadSequence("resources/yos/t.txt");
    seq = loadSequence("resources/bamboo_2gray/t.txt");
   //seq = loadSequence("resources/gsalesman/t.txt");
      CMatrix<float> img1;
      CMatrix<float> img2;
      bool presmoothing=false;
       CTensor<float> opticFlow;

    for (int i = 0; i < seq.size()-1; ++i){

	    img1 = seq(i);
	    img2 = seq(i+1);
        opticFlow(img1.xSize(),img1.ySize(),2);

      float alpha = 400;
      //float treshold= 0.0000000000001;
      float treshold= 0.0000001;
        opticFlow=Horn_SchunkOptFlow(img1 ,img2, sigma,  presmoothing, alpha,  treshold, method_choice);

    CTensor<float> Horn_SchunkFlowRGB(img1.xSize(), img1.ySize(),3);

    flowToImage(opticFlow, Horn_SchunkFlowRGB);

	Horn_SchunkFlowRGB.writeToPPM(("result/"+folderNameInput + std::to_string(i)+ ".ppm").c_str());

    }

  return 0;
}

