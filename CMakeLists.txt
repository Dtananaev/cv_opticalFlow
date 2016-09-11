cmake_minimum_required(VERSION 2.8)

project(optic_flow)

#Box filter
add_executable (Lucas_Kandade Lucas_Kandade/Lucas_Kandade.cpp )
target_link_libraries (Lucas_Kandade)

 
