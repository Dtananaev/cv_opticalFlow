/*
 * File: Brox_2004.cpp
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */

#include "optic_flow.h"
//#include "boundary.h"
#include <math.h>       /* modf */
#include "CMatrix.h"
#include "CFilter.h" // Csmooth and Cderivative
#include "CTensor.h"
#include "load_sequence.h"
#include "flowToImage.h"
#include "string"

//Cut  boundaries
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


//  diff -0,5 0 0,5
void diffXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy){

       int w =Image.xSize();
       int h =Image.ySize();
         dx.setSize(Image.xSize(),Image.ySize());
         dy.setSize(Image.xSize(),Image.ySize());
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

CTensor<float> warpImage(CMatrix<float> image1,CMatrix<float> image2, CMatrix<float> u, CMatrix<float> v){
        
    CMatrix<float> warpedImage(image1.xSize(), image1.ySize(),0);
    CMatrix<float> labels(image1.xSize(), image1.ySize(),1); //need to  make the values outside the image invalid
    CTensor<float> result(image1.xSize(), image1.ySize(),2);
   
    double  fractpart_u, intpart_u,fractpart_v, intpart_v;

  //  std::cout<<"image1 size x"<<image1.xSize()<<"\n";
  //  std::cout<<"image1 size y"<<image1.ySize()<<"\n";

   // std::cout<<"image2 size x"<<image2.xSize()<<"\n";
    //std::cout<<"image2 size y"<<image2.ySize()<<"\n";


   // std::cout<<"u size x"<<u.xSize()<<"\n";
  //  std::cout<<"u size y"<<u.ySize()<<"\n";

  //  std::cout<<"v size x"<<v.xSize()<<"\n";
  //  std::cout<<"v size y"<<v.ySize()<<"\n";
 //   std::cin.get();
    for(int x=0; x<image1.xSize();x++)
        for(int y=0; y<image1.ySize();y++){
             // split optic flow values on integer and fractal parts
            fractpart_u = modf (u(x,y) , &intpart_u);
            fractpart_v = modf (v(x,y) , &intpart_v);

             //compute the indices for diagonal pixels x1,y1 and x2,y2 of the pixels for bilinear interpolation

                /*
                        pix(x1,y1)  pix(x2,y1)

                        pix(x1,y2)  pix(x2,y2)
    
                */
            int x1,x2,y1,y2;

              x1=x+intpart_u;
              y1=y+intpart_v;
             if(fractpart_u < 0.5){
                x2=x-1 + intpart_u; 
            }  else{ x2=x+1 + intpart_u;} 
              if(fractpart_v < 0.5){
                y2=y-1 + intpart_v; 
            }  else{ y2=y+1 + intpart_v;}  

            //get the gray value from 4 pixels for interpolation from img2 if they inside image boundary otherwise use value from x,y pixel from image1 but set label matrix 0 in order to block this value for optic flow computation for image term (not for smoothsness)
            /* New notation
            A1 A2
            B1 B2
            */
            float A1,A2,B1,B2;
            //For A1 
            if( x1>=0 && x1<image2.xSize() && y1>=0 && y1<image2.ySize() ){
                A1=image2(x1,y1);
             }else{A1=image1(x,y); labels(x,y)=0;}

            //For A2
            if( x2>=0 && x2<image2.xSize() && y1>=0 && y1<image2.ySize() ){
                A2=image2(x2,y1);
             }else{A2=image1(x,y); labels(x,y)=0; }


             //For B1
            if( x1>=0 && x1<image2.xSize() && y2>=0 && y2<image2.ySize() ){
                B1=image2(x1,y2);
             }else{B1=image1(x,y); labels(x,y)=0; }   
  
              //For B2
            if( x2>=0 && x2<image2.xSize() && y2>=0 && y2<image2.ySize() ){
                B2=image2(x2,y2);
             }else{B2=image1(x,y); labels(x,y)=0; }                   
            
            
            
            //Bilinear interpolation  
             /*
                ai=(1-fractpart_u)*A1 +fractpart_u*A2
                ai_1=(1-fractpart_u)*B1 +fractpart_u*B2

                new_value=(1-fractpart_v)*ai+fractpart_v*ai_1

            */       
                
            float ai=(1-fractpart_u)*A1 +fractpart_u*A2;
            float ai_1=(1-fractpart_u)*B1 +fractpart_u*B2;

            warpedImage(x,y)=(1-fractpart_v)*ai+fractpart_v*ai_1;

    }
    
    result.putMatrix(warpedImage,0);
    result.putMatrix(labels,1);
 return result;
}


//Gauss-Seidel method
CTensor<float> Brox(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold, int lvl, float downsamplStep){

         /*Gauss-Seidel method
            Linear system:
            Ax=b;  Split A=D+L+U; D-diagonal part of matrix ; L- lower diagonal part; U - upper diagonal part
            Iterative equation:
            x(k+1)=D^-1 * (b-L*x(k+1)+U*x(k)) 
        */
    
   
    int width =image1.xSize();
    int height=image1.ySize();

    CTensor<float> result(width,height,2); //u and v optic flow result


//INITIALIZATION 

    //first downsample both images the coarse lvl
    //assume coarse lvl is 20 with step of reduction 0.95 then
    double currentLVL = pow(downsamplStep, double(lvl));
     
    int aNewXSize=trunc(width*currentLVL);
    int aNewYSize= trunc(height*currentLVL);     
        std::cout<<"init aNewXSize   "<<aNewXSize<<"\n";
        std::cout<<"init aNewYSize   "<<aNewYSize<<"\n";  
     
    CMatrix<float> CoarseImage1(image1);
    CMatrix<float> CoarseImage2(image2);

  CoarseImage1.downsampleBilinear(aNewXSize, aNewYSize);
  CoarseImage2.downsampleBilinear(aNewXSize, aNewYSize);

  // apply derivatives on downsampled images
    CMatrix<float> Ix(aNewXSize,aNewYSize);
    CMatrix<float> Iy(aNewXSize,aNewYSize); 
    diffXY(CoarseImage1, Ix, Iy);


    CMatrix<float> Iz(CoarseImage2);  
      for(int x=0; x<aNewXSize; x++)
    	for(int y=0;y<aNewYSize;y++) { 
        Iz(x,y)-=CoarseImage1(x,y);
}

    //Execute optic flow on the coarsed images for initialization
      bool presmoothing=false;
      CTensor<float>  opticFlow(aNewXSize,aNewYSize,2);
      OF flow;//init optic flow class
     // float gamma =10;
        //optic flow execut
        opticFlow=flow.Horn_SchunkOptFlow(CoarseImage1 ,CoarseImage2, 1,  presmoothing, alpha,  treshold, 10);
 
         //separate optic flow on the u and v matrices (just for convinience)
        CMatrix<float> u(opticFlow.getMatrix(0));
        CMatrix<float> v(opticFlow.getMatrix(1));
 
     //CTensor<float> Horn_SchunkFlowRGB(CoarseImage1.xSize(), CoarseImage2.ySize(),3);
   // flowToImage(opticFlow, Horn_SchunkFlowRGB);
	//Horn_SchunkFlowRGB.writeToPPM("result/Init_flow.ppm");
 // std::cin.get();
//End of initialization: we have warped image and optic flow on coarse lvl now begin iterate

    float e = 0.001; // TV regularizer
    float diff_u=10;
    float diff_v=10;
    int number_of_pixels;
 
 
   //  CMatrix<float> index; //matrix which block pixels outside of the image border for data term 
    int boundary=10;
            lvl-=1;

 for (int i=lvl; i>=0; i--) {  

          std::cout<<"lvl "<<i<<"\n";
std::cin.get();
        currentLVL = pow(downsamplStep, double(i));

        aNewXSize=trunc(width*currentLVL);
        aNewYSize= trunc(height*currentLVL);
        number_of_pixels=aNewXSize*aNewYSize;

        //upsample optic flow from previous step    
        u.upsampleBilinear(aNewXSize,aNewYSize);
        v.upsampleBilinear(aNewXSize,aNewYSize);
        CMatrix<float> du(aNewXSize,aNewYSize,0);
        CMatrix<float> dv(aNewXSize,aNewYSize,0); 
        CMatrix<float> du_new(aNewXSize,aNewYSize,0);
        CMatrix<float> dv_new(aNewXSize,aNewYSize,0);
       // CMatrix<float> u_new(aNewXSize,aNewYSize,1);
       // CMatrix<float> v_new(aNewXSize,aNewYSize,1); 
     

        //downsample images to necessary lvl
        CoarseImage1 = image1;
        CoarseImage2 = image2;

  
        CoarseImage1.downsampleBilinear(aNewXSize, aNewYSize);
        CoarseImage2.downsampleBilinear(aNewXSize, aNewYSize);
        //get warped image


        CTensor<float> warpedImage= warpImage(CoarseImage1,CoarseImage2, u, v);
        CMatrix<float> warp(warpedImage.getMatrix(0));        
        CMatrix<float>  index(warpedImage.getMatrix(1));  
        // apply derivatives on downsampled images
            diffXY(warp, Ix, Iy);
            Iz=warp;
         for(int x=0; x<aNewXSize; x++)
    	    for(int y=0;y<aNewYSize;y++) { 
             Iz(x,y)-=CoarseImage1(x,y);
          }

       //for simplicity extend borders
    Iz=Neumann_bound(Iz,boundary);
    Ix=Neumann_bound(Ix,boundary);
    Iy=Neumann_bound(Iy,boundary); 
    index=Neumann_bound(index,boundary);
    u=Neumann_bound(u,boundary);
    v=Neumann_bound(v,boundary); 

    du=Neumann_bound(du,boundary);
    dv=Neumann_bound(dv,boundary);        
    du_new=Neumann_bound(du_new,boundary);
    dv_new=Neumann_bound(dv_new,boundary);
   do{ 


      
       diff_u=0;
       diff_v=0;  



    for(int x=boundary; x<Iz.xSize()-boundary; x++)
    		for(int y=boundary;y<Iz.ySize()-boundary;y++) {  
       
        
    
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
     
        float g_1,g_2,g_3,g_4;




/*
        g_1 = 0.5* (A + B )
        g_2 = 0.5* (A + C )
        g_3 = 0.5* (A + D )   
        g_4 = 0.5* (A + E ) 
*/
   //Calculating A  = g(x,y)    

    
        int a = x;
        int b = y;

       float SQgradU =  0.25* ( u(a+1,b)-u(a-1,b) ) * ( u(a+1,b)-u(a-1,b) )  + 0.25* ( u(a,b+1)-u(a,b-1) ) * ( u(a,b+1)-u(a,b-1) );

       float SQgradV =  0.25* ( v(a+1,b)-v(a-1,b) ) * ( v(a+1,b)-v(a-1,b) )  + 0.25* ( v(a,b+1)-v(a,b-1) ) * ( v(a,b+1)-v(a,b-1) );

       float A = 1/sqrt(   SQgradU + SQgradV +e*e  );


        //Calculating B = g(x+1,y)  
         a = x+1;  b = y;

           SQgradU =  0.25* ( u(a+1,b)-u(a-1,b) ) * ( u(a+1,b)-u(a-1,b) )  + 0.25* ( u(a,b+1)-u(a,b-1) ) * ( u(a,b+1)-u(a,b-1) );

          SQgradV =  0.25* ( v(a+1,b)-v(a-1,b) ) * ( v(a+1,b)-v(a-1,b) )  + 0.25* ( v(a,b+1)-v(a,b-1) ) * ( v(a,b+1)-v(a,b-1) );
        
       float B = 1/sqrt(   SQgradU + SQgradV +e*e  );

        //Calculating C = g(x-1,y)  
         a =x-1; b = y;
 
             SQgradU =  0.25* ( u(a+1,b)-u(a-1,b) ) * ( u(a+1,b)-u(a-1,b) )  + 0.25* ( u(a,b+1)-u(a,b-1) ) * ( u(a,b+1)-u(a,b-1) );

          SQgradV =  0.25* ( v(a+1,b)-v(a-1,b) ) * ( v(a+1,b)-v(a-1,b) )  + 0.25* ( v(a,b+1)-v(a,b-1) ) * ( v(a,b+1)-v(a,b-1) );
        
       float C = 1/sqrt(   SQgradU + SQgradV +e*e  );

        //Calculating D = g(x,y+1)  
         a =x; b = y+1;
    
            SQgradU =  0.25* ( u(a+1,b)-u(a-1,b) ) * ( u(a+1,b)-u(a-1,b) )  + 0.25* ( u(a,b+1)-u(a,b-1) ) * ( u(a,b+1)-u(a,b-1) );
 
          SQgradV =  0.25* ( v(a+1,b)-v(a-1,b) ) * ( v(a+1,b)-v(a-1,b) )  + 0.25* ( v(a,b+1)-v(a,b-1) ) * ( v(a,b+1)-v(a,b-1) );

       float D = 1/sqrt(   SQgradU + SQgradV +e*e  );

        //Calculating E = g(x,y-1)  
         a =x; b =y-1;
 
             SQgradU =  0.25* ( u(a+1,b)-u(a-1,b) ) * ( u(a+1,b)-u(a-1,b) )  + 0.25* ( u(a,b+1)-u(a,b-1) ) * ( u(a,b+1)-u(a,b-1) );

          SQgradV =  0.25* ( v(a+1,b)-v(a-1,b) ) * ( v(a+1,b)-v(a-1,b) )  + 0.25* ( v(a,b+1)-v(a,b-1) ) * ( v(a,b+1)-v(a,b-1) );
        
       float E = 1/sqrt(   SQgradU + SQgradV +e*e  );


    //Calculating our functions for solver and keep them constant for one iteration
       g_1 = 0.5* (A + B );
       g_2 = 0.5* (A + C );
       g_3 = 0.5* (A + D );   
       g_4 = 0.5* (A + E );

    //calculate constant nonlinear part
   float Ct=1/sqrt(Iz(x,y)*Iz(x,y)+e*e);  
   //Gauss-Seidel linear solver for Optic flow
   float SmoothConst=g_1+g_2+g_3+g_4;
  float SmoothU= g_1*u(x+1,y)+g_2*u(x-1,y)+ g_3*u(x,y+1)+g_4*u(x,y-1);
   float SmoothV= g_1*v(x+1,y)+g_2*v(x-1,y)+ g_3*v(x,y+1)+g_4*v(x,y-1);
   //float SmoothU= g_1*u(x+1,y)+g_2*u(x-1,y)+ g_3*u(x,y+1)+g_4*u(x,y-1) +g_1*du(x+1,y)+g_2*du(x-1,y)+ g_3*du(x,y+1)+g_4*du(x,y-1);
  // float SmoothV= g_1*v(x+1,y)+g_2*v(x-1,y)+ g_3*v(x,y+1)+g_4*v(x,y-1) + g_1*v(x+1,y)+g_2*dv(x-1,y)+ g_3*dv(x,y+1)+g_4*dv(x,y-1);
   float denomU= Ct*Ix(x,y)*Ix(x,y)*index(x,y)+alpha*SmoothConst;
   float denomV= Ct*Iy(x,y)*Iy(x,y)*index(x,y)+alpha*SmoothConst;
      

     du_new(x,y)= ( alpha*SmoothU - SmoothConst* u(x,y) - ( Ct*Iy(x,y)*Ix(x,y)*dv(x,y)+Iz(x,y)*Ix(x,y) ) * index(x,y) )/denomU;


      dv_new(x,y)= ( alpha*SmoothV - SmoothConst* v(x,y) - ( Ct*Iy(x,y)*Ix(x,y)*du(x,y)+Iz(x,y)*Iy(x,y) ) * index(x,y) )/denomV;  
   
   //du_new(x,y)= (1 / ((g_1+g_2+g_3+g_4)*alpha + Ix(x,y)*Ix(x,y)) ) * ( alpha* (g_2*du(x-1,y)+g_1*du(x+1,y)+g_4*du(x,y-1)+g_3*du(x,y+1))  -Ix(x,y)*Iy(x,y)*dv(x,y) - Ix(x,y)*Iz(x,y) );
            

       //     dv_new(x,y) = (1/ ((g_1+g_2+g_3+g_4)*alpha+ Iy(x,y)*Iy(x,y)) )*(alpha* (g_2*dv(x-1,y)+g_1*dv(x+1,y)+g_4*dv(x,y-1)+g_3*dv(x,y+1)) - Iy(x,y)*Ix(x,y)*du(x,y) -Iy(x,y)*Iz(x,y) );
          

        diff_u+=  (du_new(x,y)- du(x,y))* (du_new(x,y)- du(x,y));
        diff_v+=  (dv_new(x,y)- dv(x,y))* (dv_new(x,y)- dv(x,y));    

     
    }
 
        du=du_new;
        dv=dv_new;
        diff_u=  diff_u/ number_of_pixels;
        diff_v=  diff_v/ number_of_pixels;
        std::cout<<"diff_u "<<diff_u<<"\n";
        std::cout<<"diff_v "<<diff_v<<"\n";

    }  while ( diff_u > treshold && diff_v> treshold);

    for(int x=boundary; x<Iz.xSize()-boundary; x++)
    		for(int y=boundary;y<Iz.ySize()-boundary;y++) {  
        u(x,y)=u(x,y)+du(x,y);
        v(x,y)=v(x,y)+dv(x,y);

    }
        u=cut(u,boundary);
        v=cut(v,boundary);  
       
}


       result.putMatrix(u,0);      
       result.putMatrix(v,1);      

    return result; 

}



int main(int argc, char** argv) {  

  

    //optic flow parameters
      float alpha = 500;
     // float gamma =10;
      float treshold= 0.00000001;
  //  float treshold= -0.00000001;
        int lvl=20; 
        float downsamplStep=0.95;
    std::string folderNameInput;

	//int sigma=1;



	if (argc==2){
		folderNameInput=argv[1];

	//} else if (argc==3){
	//	folderNameInput=argv[1];
		//sigma=atoi(argv[2]);

	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		std::cout<<"Usage: Lucas_Kanade inputfolder <variance of Gauss kernell>"<<"\n";
		std::cout<<"The command should contain at least input file name. The default Gauss kernel sigma=1."<<"\n";
		return 0;    
	}
    

    CVector<CMatrix<float> > seq;
    
   // seq = loadSequence("resources/cropped-street/t.txt");
   seq = loadSequence("resources/yos/t.txt");
   //seq = loadSequence("resources/gsalesman/t.txt");
      CMatrix<float> img1;
      CMatrix<float> img2;
       CTensor<float> opticFlow;
    for (int i = 0; i < seq.size()-1; ++i){
        //get images from sequence
	    img1 = seq(i);
	    img2 = seq(i+1);
        opticFlow(img1.xSize(),img1.ySize(),2);


    opticFlow= Brox(img1,  img2, alpha,  treshold,  lvl,  downsamplStep);

     
        //start calculateincrement  of the optic flow by the Brox method


    //visualize resulted warped image       
    //CMatrix<float> warpedImageVis(warpedImage.getMatrix(0));  
   // warpedImageVis.writeToPGM(("result/"+folderNameInput + std::to_string(i)+ "_warp.pgm").c_str());

    //visualize optic flow 
    CTensor<float> Horn_SchunkFlowRGB(img1.xSize(), img1.ySize(),3);
    flowToImage(opticFlow, Horn_SchunkFlowRGB);
	Horn_SchunkFlowRGB.writeToPPM(("result/"+folderNameInput + std::to_string(i)+ ".ppm").c_str());

    }



    
  return 0;
}
