/*
 * File: Brox_2004.cpp
 * 
 * Author: Denis Tananaev
 *
 * Date: 18.09.2016
 */

#include "optic_flow.h"



int main(int argc, char** argv) {  

    OF flow;
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
    

    CVector<CMatrix<float> > seq;
    
   // seq = loadSequence("resources/cropped-street/t.txt");
   seq = loadSequence("resources/yos/t.txt");
   //seq = loadSequence("resources/gsalesman/t.txt");
      CMatrix<float> img1;
      CMatrix<float> img2;
      bool presmoothing=true;
       CTensor<float> opticFlow;

    for (int i = 0; i < seq.size()-1; ++i){

	    img1 = seq(i);
	    img2 = seq(i+1);
        opticFlow(img1.xSize(),img1.ySize(),2);

      float alpha = 500;
      float gamma =20;
      float treshold= 0.00000001;
 
        opticFlow=flow.Horn_SchunkOptFlow(img1 ,img2, sigma,  presmoothing, alpha,  treshold, gamma);

    CTensor<float> Horn_SchunkFlowRGB(img1.xSize(), img1.ySize(),3);

    flowToImage(opticFlow, Horn_SchunkFlowRGB);

	Horn_SchunkFlowRGB.writeToPPM(("result/"+folderNameInput + std::to_string(i)+ ".ppm").c_str());

    }

  return 0;
}
