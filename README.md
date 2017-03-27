Computer vision: Optical Flow
====================================================
All Optical Flow implementations can process only files in pgm format.
In order to visualize the optical flow the next color scheme is used:
<p align="center">
  <img src="https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/Color_scheme.JPG" width="700"/>
</p>


[![Build Status](https://travis-ci.org/Dtananaev/cv_opticalFlow.svg?branch=master)](https://travis-ci.org/Dtananaev/cv_opticalFlow)
[![BSD2 License](http://img.shields.io/badge/license-BSD2-brightgreen.svg)](https://github.com/Dtananaev/cv_opticFlow/blob/master/LICENSE.md) 
     
It contains:

* Lucas_Kanade - Lucas-Kanade optical flow. One of the first successful implementation of the optical flow based on gray value constancy assumption. 
 [![Lucas_Kanade](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/Lucas_Kanade.JPG)](https://www.youtube.com/watch?v=wd3EbR8unJQ)
     * To install use in terminal: 
       * cd ../cv_opticFlow/Lucas_Kanade
       * make
     * To run: ./Lucas_Kanade name_of_folder_with_picture_sequence 
* Horn_Schunck - Horn_Schunck optical flow based on variational methods.
 [![Horn_Schunck](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/Horn_Schunck.JPG)](https://www.youtube.com/watch?v=vQioi02NS9A)
     * To install use in terminal: 
       * cd ../cv_opticFlow/Horn_Schunck
       * make
     * To run: ./Horn_Schunck name_of_folder_with_picture_sequence 
     
* Horn_Schunck_TV - Horn_Schunck optical flow based on variational methods with total variation (TV) smoothness term (robust penalizer for making sharp edges).
 [![ Horn_Schunck_TV ](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/HS_TV.JPG)](https://www.youtube.com/watch?v=hO7HGA_PFD8)
     * To install use in terminal: 
       * cd ../cv_opticFlow/Horn_Schunck_TV
       * make
     * To run: ./Horn_Schunck_TV name_of_folder_with_picture_sequence 
     
* Horn_Schunck_TVbothTerms - Horn_Schunck optical flow based on variational methods with total variation (TV) both terms.
 [![ Horn_Schunck_TVbothTerms ](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/HS_both.JPG)](https://www.youtube.com/watch?v=bse2mM_eRr4)
     * To install use in terminal: 
       * cd ../cv_opticFlow/Horn_Schunck_TVbothTerms
       * make
     * To run: ./Horn_Schunck_TVboth name_of_folder_with_picture_sequence 
     
* Horn_Schunck_TVbothTerms_Grad - Horn_Schunck optical flow based on variational methods with total variation (TV) both terms  plus gradient constancy assumption which makes optical flow robust to illumination changes.
 [![ Horn_Schunck_TVbothTerms_Grad ](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/HS_grad.JPG)](https://www.youtube.com/watch?v=eNuUIlOj4SA)
     * To install use in terminal: 
       * cd ../cv_opticFlow/Horn_Schunck_TVbothTerms_Grad
       * make
     * To run: ./Horn_Schunck_TVbothTerms_Grad name_of_folder_with_picture_sequence 
     
* Brox_2004 - Thomas Brox large displacements optical flow  based on variational methods without gradient constancy assumption. 
 [![ Brox_2004 ](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/Brox.JPG)](https://www.youtube.com/watch?v=I7h4qt0taXA&list=PLm9UBVQa6prnXGdCuKKLre8-yAyfAPHi-&index=1)
     * To install use in terminal: 
       * cd ../cv_opticFlow/Brox_2004
       * make
     * To run: ./Brox_2004 name_of_folder_with_picture_sequence 
     
* denoiz - filtering algorithm equivalent to Gauss blur filtering but based on the energy minimization approach. Could be helpfull for improving quality of the optical flow for noisy images.
      * To install use in terminal: 
       * cd ../cv_opticFlow/denoiz
       * make
     * To run: ./denoiz name_of_file
 
