/*
 * File: optic_flow.cpp
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */




#include "optic_flow.h"


OF::OF(){
}

OF::~OF(){
}
//Dirichlet boundary conditions
CMatrix<float> OF::Dirichlet_bound(CMatrix<float> aImage,int border_size){

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
//Cut  boundaries
CMatrix<float> OF::cut(CMatrix<float>& image,int border_size){ 

	CMatrix<float> realimage(image.xSize()-2*border_size,image.ySize()-2*border_size);
	for(int x=0;x<realimage.xSize();x++)
		for(int y=0; y<realimage.ySize();y++)
		{
			realimage(x,y)=image(x+border_size,y+border_size);
		}

		return realimage;
}

// Neumann boundry conditions
CMatrix<float> OF::Neumann_bound(CMatrix<float> aImage,int border_size){

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
CMatrix<float> OF::Gauss(int sigma){

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
CMatrix<float> OF::Gfilter(CMatrix<float> Gauss, CMatrix<float> boundary_Image,int border){
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
void OF::diffXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy){
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
void OF::difForwXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy){
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

void OF::calculateTV(CMatrix<float> u_k, CMatrix<float> v_k, CMatrix<float>& g_1,CMatrix<float>& g_2,CMatrix<float>& g_3,CMatrix<float>& g_4, int border ){
     float e = 0.0001; // TV regularizer

 //Lagged diffusivity scheme: First calculate the nonlinear part which is TV diffusion
         /*    

        //Calculating the constant functions which represents the gradient flow between two pixels 
         g_1 = g(x+0,5,y) = 0.5 * (g(x+1,y)+g(x,y))
         g_2 = g(x-0,5,y) = 0.5 * (g(x-1,y)+g(x,y))
         g_3 = g(x,y+0,5) = 0.5 * (g(x,y+1)+g(x,y))
         g_4 = g(x,y-0,5) = 0.5 * (g(x,y-1)+g(x,y))

        //Here g functions is represented by the TV for u and v optic flow vectors (smoothsness term of the Horn-Schunck optic flow). The TV function:
        //Lets make new notation
       A =  g(x,y) =  1/ sqrt( grad(u(x,y))^2   + grad(v(x,y))^2 + e^2)
       B = g(x+1,y) =  1/ sqrt( grad(u(x+1,y))^2   + grad(v(x+1,y))^2 + e^2)
       C = g(x-1,y) =  1/ sqrt( grad(u(x-1,y))^2   + grad(v(x-1,y))^2 + e^2)
       D = g(x,y+1) = 1/ sqrt( grad(u(x,y+1))^2   + grad(v(x,y+1))^2 + e^2)
       E = g(x,y-1) = 1/ sqrt( grad(u(x,y-1))^2   + grad(v(x,y-1))^2 + e^2)

        from here 

        g_1 = 0.5* (A + B )
        g_2 = 0.5* (A + C )
        g_3 = 0.5* (A + D )   
        g_4 = 0.5* (A + E ) 
  
        */
 for(int x=border; x< u_k.xSize()-border; x++)
    		for(int y=border;y< u_k.ySize()-border;y++) { 


        int a = x;
        int b = y;

       float SQgradU =  0.25* ( u_k(a+1,b)-u_k(a-1,b) ) * ( u_k(a+1,b)-u_k(a-1,b) )  + 0.25* ( u_k(a,b+1)-u_k(a,b-1) ) * ( u_k(a,b+1)-u_k(a,b-1) );

       float SQgradV =  0.25* ( v_k(a+1,b)-v_k(a-1,b) ) * ( v_k(a+1,b)-v_k(a-1,b) )  + 0.25* ( v_k(a,b+1)-v_k(a,b-1) ) * ( v_k(a,b+1)-v_k(a,b-1) );

       float A = 1/sqrt(   SQgradU + SQgradV +e*e  );


        //Calculating B = g(x+1,y)  
         a = x+1;  b = y;

        SQgradU =  0.25* ( u_k(a+1,b)-u_k(a-1,b) ) * ( u_k(a+1,b)-u_k(a-1,b) )  + 0.25* ( u_k(a,b+1)-u_k(a,b-1) ) * ( u_k(a,b+1)-u_k(a,b-1) );

        SQgradV =  0.25* ( v_k(a+1,b)-v_k(a-1,b) ) * ( v_k(a+1,b)-v_k(a-1,b) )  + 0.25* ( v_k(a,b+1)-v_k(a,b-1) ) * ( v_k(a,b+1)-v_k(a,b-1) );
        
       float B = 1/sqrt(   SQgradU + SQgradV +e*e  );

        //Calculating C = g(x-1,y)  
         a =x-1; b = y;
 
        SQgradU =  0.25* ( u_k(a+1,b)-u_k(a-1,b) ) * ( u_k(a+1,b)-u_k(a-1,b) )  + 0.25* ( u_k(a,b+1)-u_k(a,b-1) ) * ( u_k(a,b+1)-u_k(a,b-1) );

        SQgradV =  0.25* ( v_k(a+1,b)-v_k(a-1,b) ) * ( v_k(a+1,b)-v_k(a-1,b) )  + 0.25* ( v_k(a,b+1)-v_k(a,b-1) ) * ( v_k(a,b+1)-v_k(a,b-1) );
        
       float C = 1/sqrt(   SQgradU + SQgradV +e*e  );

        //Calculating D = g(x,y+1)  
         a =x; b = y+1;
    
        SQgradU =  0.25* ( u_k(a+1,b)-u_k(a-1,b) ) * ( u_k(a+1,b)-u_k(a-1,b) )  + 0.25* ( u_k(a,b+1)-u_k(a,b-1) ) * ( u_k(a,b+1)-u_k(a,b-1) );

        SQgradV =  0.25* ( v_k(a+1,b)-v_k(a-1,b) ) * ( v_k(a+1,b)-v_k(a-1,b) )  + 0.25* ( v_k(a,b+1)-v_k(a,b-1) ) * ( v_k(a,b+1)-v_k(a,b-1) );
        
       float D = 1/sqrt(   SQgradU + SQgradV +e*e  );

        //Calculating E = g(x,y-1)  
         a =x; b =y-1;
 
        SQgradU =  0.25* ( u_k(a+1,b)-u_k(a-1,b) ) * ( u_k(a+1,b)-u_k(a-1,b) )  + 0.25* ( u_k(a,b+1)-u_k(a,b-1) ) * ( u_k(a,b+1)-u_k(a,b-1) );

        SQgradV =  0.25* ( v_k(a+1,b)-v_k(a-1,b) ) * ( v_k(a+1,b)-v_k(a-1,b) )  + 0.25* ( v_k(a,b+1)-v_k(a,b-1) ) * ( v_k(a,b+1)-v_k(a,b-1) );
        
       float E = 1/sqrt(   SQgradU + SQgradV +e*e  );


    //Calculating our functions for solver and keep them constant for one iteration
       g_1(x,y) = 0.5* (A + B );
       g_2(x,y) = 0.5* (A + C );
       g_3(x,y) = 0.5* (A + D );   
       g_4(x,y) = 0.5* (A + E );
    
    }


}


//Gauss-Seidel method
CTensor<float> OF::GaussSeidelHSTV(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold, float gamma){
    int counter=0;
  float e = 0.0001; // TV regularizer

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

//apply second derivative
    CMatrix<float> Ixx(width,height);
    CMatrix<float> Iyy(width,height);    
    CMatrix<float> Ixy(width,height);
    CMatrix<float> Iyz(width,height); 
    CMatrix<float> Ixz(width,height); 
    diffXY(Ix, Ixx, Ixy);
    diffXY(Iy, Ixy, Iyy);
    diffXY(Iz, Ixz, Iyz);

    //To make computation simpler  and have all matrices of the same size apply Dirichlet boundary conditions with border size 1 to derivatives and image1
    int border=3;
    image1=Neumann_bound(image1,border);
    Iz=Neumann_bound(Iz,border);
    Ix=Neumann_bound(Ix,border);
    Iy=Neumann_bound(Iy,border); 
  
    Ixz=Neumann_bound(Ixz,border);
    Iyz=Neumann_bound(Iyz,border);
    Ixx=Neumann_bound(Ixx,border);
    Ixy=Neumann_bound(Ixy,border);  
    Iyy=Neumann_bound(Iyy,border);

     CMatrix<float> g_1(image1.xSize(),image1.ySize(),1);
     CMatrix<float> g_2(image1.xSize(),image1.ySize(),1);
     CMatrix<float> g_3(image1.xSize(),image1.ySize(),1);
     CMatrix<float> g_4(image1.xSize(),image1.ySize(),1);
    //Horn-Schunck optic flow
     //CMatrix<float> u_k(image1);
    // CMatrix<float> v_k(image1);     
    CMatrix<float> u_k(image1.xSize(),image1.ySize(),100);
     CMatrix<float> v_k(image1.xSize(),image1.ySize(),100);
     CMatrix<float> u_k_new(image1.xSize(),image1.ySize(),0);
     CMatrix<float> v_k_new(image1.xSize(),image1.ySize(),0);
     float diff_u,diff_v;

   
                    
  do{ 
       diff_u=0;
       diff_v=0;  
    for(int x=border; x<image1.xSize()-border; x++)
    		for(int y=border;y<image1.ySize()-border;y++) {  
        

// Now calculate the gray value cosntacy assumption nonlinearity

        float C=1/sqrt( (Ix(x,y)*u_k(x,y) + Iy(x,y)*v_k(x,y) + Iz(x,y))*(Ix(x,y)*u_k(x,y) + Iy(x,y)*v_k(x,y) + Iz(x,y)) +e*e );

//Now calculate the gradient value constancy assumption

      float gradC=gamma*1/sqrt(  pow( (Ixx(x,y)*u_k(x,y) + Ixy(x,y)* v_k(x,y) + Ixz(x,y)),2) +  pow( (Ixy(x,y)*u_k(x,y) + Iyy(x,y)* v_k(x,y) + Iyz(x,y)),2) + e*e  );

            //Gauss-Seidel linear solver


            u_k_new(x,y)= (1 / ((g_1(x,y)+g_2(x,y)+g_3(x,y)+g_4(x,y))*alpha + C*Ix(x,y)*Ix(x,y) + gradC*(Ixx(x,y)*Ixx(x,y)+Ixy(x,y)*Ixy(x,y)) ) ) * ( alpha* (g_2(x,y)*u_k_new(x-1,y)+g_1(x,y)*u_k(x+1,y)+g_4(x,y)*u_k_new(x,y-1)+g_3(x,y)*u_k(x,y+1))  -C*Ix(x,y)*Iy(x,y)*v_k(x,y) - C*Ix(x,y)*Iz(x,y) - gradC*v_k(x,y)*(Ixx(x,y)*Ixy(x,y)+ Iyy(x,y)*Ixy(x,y)) -gradC*(Ixz(x,y)*Ixx(x,y)+Iyz(x,y)*Ixy(x,y) )  );
            

            v_k_new(x,y) = (1/ ((g_1(x,y)+g_2(x,y)+g_3(x,y)+g_4(x,y))*alpha+ C*Iy(x,y)*Iy(x,y) + gradC*(Iyy(x,y)*Iyy(x,y)+Ixy(x,y)*Ixy(x,y))    ) )*(alpha* (g_2(x,y)*v_k_new(x-1,y)+g_1(x,y)*v_k(x+1,y)+g_4(x,y)*v_k_new(x,y-1)+g_3(x,y)*v_k(x,y+1)) - C*Iy(x,y)*Ix(x,y)*u_k_new(x,y) -C*Iy(x,y)*Iz(x,y)    - gradC*u_k_new(x,y)*(Ixx(x,y)*Ixy(x,y)+ Iyy(x,y)*Ixy(x,y)) -gradC*(Ixz(x,y)*Ixy(x,y)+Iyz(x,y)*Iyy(x,y)  ) );

        diff_u+= (u_k_new(x,y)-u_k(x,y))*(u_k_new(x,y)-u_k(x,y));
        diff_v+= (v_k_new(x,y)-v_k(x,y))*(v_k_new(x,y)-v_k(x,y));    
            


    }
         u_k=u_k_new;
         v_k=v_k_new;  
        // u_k=cut(u_k,border);
        // v_k=cut(v_k,border); 
         //u_k=Neumann_bound(u_k,border);
        // v_k=Neumann_bound(v_k,border);


       // calculateTV(u_k,v_k, g_1,g_2, g_3, g_4, border );
          diff_u=  diff_u/ number_of_pixels;
          diff_v=  diff_v/ number_of_pixels;
        std::cout<<"diff_u "<<diff_u<<"\n";
         std::cout<<"diff_v "<<diff_v<<"\n";
        counter+=1;
         std::cout<<"counter "<<counter<<"\n";
        if(counter>=5000){break;}
} while ( diff_u >treshold && diff_v> treshold );

        u_k_new=cut(u_k_new,border);
        v_k_new=cut(v_k_new,border); 

       result.putMatrix(u_k_new,0);      
       result.putMatrix(v_k_new,1);      

    return result; 
}


CTensor<float> OF::Horn_SchunkOptFlow(CMatrix<float> image1, CMatrix<float> image2, int sigma, bool presmoothing, float alpha, float treshold, float gamma){
    
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
        result=GaussSeidelHSTV( image1,  image2,  alpha,  treshold,gamma);

return result;
  
}




