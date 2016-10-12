Computer vision: Optical Flow
====================================================
All Optical Flow implementations can process only files in pgm format.

[![Build Status](https://travis-ci.org/Dtananaev/cv_filters.svg?branch=master)](https://travis-ci.org/Dtananaev/cv_filters)
[![BSD2 License](http://img.shields.io/badge/license-BSD2-brightgreen.svg)](https://github.com/Dtananaev/cv_filters/blob/master/LICENSE.md) 
     
It contains:

* Lucas_Kanade - Lucas-Kanade optical flow. One of the first successful implementation of the optical flow based on gray value constancy assumption. 
 [![Lucas_Kanade](https://github.com/Dtananaev/cv_opticFlow/blob/master/pictures/Lucas_Kanade.JPG)](https://www.youtube.com/watch?v=wd3EbR8unJQ)
     * To install use in terminal: 
          * cd ../cv_opticFlow/Lucas_Kanade
          * make
     * To run: ./Lucas_Kanade name_of_folder_with_picture_sequence 
