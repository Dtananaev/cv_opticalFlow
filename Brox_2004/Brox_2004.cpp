/*
 * File: Brox_2004.cpp
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */

#include "optic_flow.h"
#include "boundary.h"
#include <math.h>       /* modf */
#include "CMatrix.h"
#include "CFilter.h" // Csmooth and Cderivative
#include "CTensor.h"
#include "load_sequence.h"
#include "flowToImage.h"
#include "string"
//  diff -0,5 0 0,5
void diffXY(CMatrix<float> Image,  CMatrix<float> &dx, CMatrix<float> &dy){
       boundary bound;
       int w =Image.xSize();
       int h =Image.ySize();
       Image=bound.Neumann_bound(Image,1);
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
             }else{A1=image1(x,y); labels(x,y)=0; }

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
    result.putMatrix(warpedImage,1);
 return result;
}


//Gauss-Seidel method
CTensor<float> Brox(CMatrix<float> image1, CMatrix<float> image2, float alpha, float treshold, int lvl, float downsamplStep){
     boundary bound;  //boundary conditions    
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
      CTensor<float> opticFlow;
      CTensor<float> warpedImage;
      OF flow;//init optic flow class
     // float gamma =10;
       
        opticFlow(aNewXSize,aNewYSize,2);


        //optic flow execut
        opticFlow=flow.Horn_SchunkOptFlow(CoarseImage1 ,CoarseImage2, 1,  presmoothing, alpha,  treshold, 10);
        std::cout<<"1"<<"\n";
         //separate optic flow on the u and v matrices (just for convinience)
        CMatrix<float> u(opticFlow.getMatrix(0));
        CMatrix<float> v(opticFlow.getMatrix(1));
 
     std::cout<<"2"<<"\n";

//End of initialization: we have warped image and optic flow on coarse lvl now begin iterate

    float e = 0.001; // TV regularizer
    float diff_u,diff_v;
    int number_of_pixels;
     CMatrix<float> du;
     CMatrix<float> dv;  
     CMatrix<float> u_new;
     CMatrix<float> v_new;  
     CMatrix<float> index; //matrix which block pixels outside of the image border for data term 
    int boundary=3;
     
  do{ 
       diff_u=0;
       diff_v=0;  
        u=bound.cut(u,boundary);
        v=bound.cut(v,boundary);
       //upsample the optic flow for one lvl
        lvl-=1;
        currentLVL = pow(downsamplStep, double(lvl));
     
        aNewXSize=trunc(width*currentLVL);
        aNewYSize= trunc(height*currentLVL);
        number_of_pixels=aNewXSize*aNewYSize;
        //upsample optic flow from previous step    
        u.upsampleBilinear(aNewXSize,aNewYSize);
        v.upsampleBilinear(aNewXSize,aNewYSize);

        //downsample images to necessary lvl
        CoarseImage1 = image1;
        CoarseImage2 = image2;

        CoarseImage1.downsampleBilinear(aNewXSize, aNewYSize);
        CoarseImage2.downsampleBilinear(aNewXSize, aNewYSize);

        //get warped image
        warpedImage(aNewXSize,aNewYSize,2);
        warpedImage= warpImage(CoarseImage1,CoarseImage2, u, v);

        Iz=warpedImage.getMatrix(0);  
        index=warpedImage.getMatrix(1);  
     // apply derivatives on downsampled images
         Ix(aNewXSize,aNewYSize);
         Iy(aNewXSize,aNewYSize); 
        diffXY(CoarseImage1, Ix, Iy);


         for(int x=0; x<aNewXSize; x++)
    	    for(int y=0;y<aNewYSize;y++) { 
             Iz(x,y)-=CoarseImage1(x,y);
          }


    //for simplicity extend borders


    Iz=bound.Neumann_bound(Iz,boundary);
    Ix=bound.Neumann_bound(Ix,boundary);
    Iy=bound.Neumann_bound(Iy,boundary); 
    u=bound.Neumann_bound(u,boundary);
    v=bound.Neumann_bound(v,boundary); 
            //reset increment
             du(Iz.xSize(),Iz.ySize());
             dv(Iz.xSize(),Iz.ySize());

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
        
       float Cd = 1/sqrt(   SQgradU + SQgradV +e*e  );

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
       g_2 = 0.5* (A + Cd );
       g_3 = 0.5* (A + D );   
       g_4 = 0.5* (A + E );

    //calculate constant nonlinear part
   float C=1/sqrt(Iz(x,y)*Iz(x,y)+e*e);
    

    // Calculate 
    
   
    //Gauss-Seidel linear solver for Optic flow
    
    

        du(x,y)=( 1/(C*Ix(x,y)*Ix(x,y) - 0.5*alpha * (g_1-g_2+g_3-g_4)) )* ( 0.5*alpha*(0.5*g_1* (u(x+2,y)-u(x,y)+u(x+1,y+1)-u(x+1,y-1) ) -0.5*g_2*(u(x,y)-u(x-2,y)+u(x-1,y+1)-u(x-1,y-1)) +0.5*g_3*(u(x+1,y+1)-u(x-1,y+1)+u(x,y+2)-u(x,y)) -0.5*g_4*(u(x+1,y-1)-u(x-1,y-1)+u(x,y)-u(x,y-2))) -( 0.5*C*Iy(x,y)*Ix(x,y)*(v(x+1,y)-v(x-1,y)+v(x,y+1)-v(x,y-1)) + C*Iz(x,y)*Ix(x,y))*index(x,y));


        dv(x,y)=( 1/(C*Iy(x,y)*Iy(x,y) - 0.5*alpha * (g_1-g_2+g_3-g_4)) )* ( 0.5*alpha*(0.5*g_1* (v(x+2,y)-v(x,y)+v(x+1,y+1)-v(x+1,y-1) ) -0.5*g_2*(v(x,y)-v(x-2,y)+v(x-1,y+1)-v(x-1,y-1)) +0.5*g_3*(v(x+1,y+1)-v(x-1,y+1)+v(x,y+2)-v(x,y)) -0.5*g_4*(v(x+1,y-1)-v(x-1,y-1)+v(x,y)-v(x,y-2))) -( 0.5*C*Iy(x,y)*Ix(x,y)*(u(x+1,y)-u(x-1,y)+u(x,y+1)-u(x,y-1)) + C*Iz(x,y)*Iy(x,y))*index(x,y));

        u_new(x,y)=u(x,y)+du(x,y);
        v_new(x,y)=v(x,y)+dv(x,y);
        
        diff_u+=  du(x,y)* du(x,y);
        diff_v+=  dv(x,y)* dv(x,y);    
            


    }
 
         u=u_new;
         v=v_new;  
      
    
          diff_u=  diff_u/ number_of_pixels;
          diff_v=  diff_v/ number_of_pixels;
        std::cout<<"diff_u "<<diff_u<<"\n";
         std::cout<<"diff_v "<<diff_v<<"\n";


} while ( diff_u >treshold && diff_v> treshold && lvl>0);

        u_new=bound.cut(u_new,boundary);
        v_new=bound.cut(v_new,boundary); 

       result.putMatrix(u_new,0);      
       result.putMatrix(v_new,1);      

    return result; 

}



int main(int argc, char** argv) {  

  

    //optic flow parameters
      float alpha = 500;
     // float gamma =10;
      float treshold= 0.00000001;
        int lvl=5; 
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
