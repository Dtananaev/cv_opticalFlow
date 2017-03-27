cmake_minimum_required(VERSION 2.8)

project(optic_flow)

include_directories(
  Lucas_Kanade/include
  Brox_2004/include  
  Horn_Schunck/include
  Horn_Schunck_TV/include  
  Horn_Schunck_TVbothTerms/include  
  Horn_Schunck_TVbothTerms_Grad/include  
)

add_executable (Lucas_Kanade Lucas_Kanade/Lucas_Kanade.cpp )
target_link_libraries (Lucas_Kanade)

add_executable (Brox_2004 Brox_2004/Brox_2004.cpp )
target_link_libraries (Brox_2004)
 
add_executable (Horn_Schunck Horn_Schunck/Horn_Schunck.cpp )
target_link_libraries (Horn_Schunck)
 
add_executable (Horn_Schunck_TV Horn_Schunck_TV/Horn_Schunck_TV.cpp )
target_link_libraries (Horn_Schunck_TV)

add_executable (Horn_Schunck_TVboth Horn_Schunck_TVbothTerms/Horn_Schunck_TVboth.cpp )
target_link_libraries (Horn_Schunck_TVboth)

add_executable (Horn_Schunck_TVbothTerms_Grad Horn_Schunck_TVbothTerms_Grad/Horn_Schunck_TVbothTerms_Grad.cpp )
target_link_libraries (Horn_Schunck_TVbothTerms_Grad)

add_executable (denoiz denoiz/denoiz.cpp )
target_link_libraries (denoiz)
