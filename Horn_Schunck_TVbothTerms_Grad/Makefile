##
#
# Author: Denis Tananaev, 2015
#
##

## Where to look for includes (default is 'here')
INCLUDE_DIRS += -I. 
INCLUDE_DIRS += -I./src
INCLUDE_DIRS += -I./include
## Compiler
CXX ?= gcc

OBJDIR=obj

## Compiler flags; extended in 'debug'/'release' rules
CXXFLAGS += -W -Wall -Wextra -Wpedantic -std=c++11 -pthread -fopenmp

## Linker flags
LDFLAGS = -pthread -lX11 -fopenmp

## Default name for the built executable
TARGET = Horn_Schunck_TVboth_grad

## Every *.c/*.cc/*.cpp file is a source file
SOURCES = $(wildcard *.c *.cc *.cpp)
## Every *.h/*.hh/*.hpp file is a header file
HEADERS = $(wildcard *.h *.hh *.hpp)

## Build a *.o object file for every source file
#OBJECTS = $(addsuffix  *.o, $(basename $(SOURCES)))
OBJECTS = $(addprefix  $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))

## Tell the 'make' program that e.g. 'make clean' is not supposed to 
#  create a file 'clean'
# 
#  "Why is it called 'phony'?" -- because it's not a real target. That 
#  is, the target name isn't a file that is produced by the commands 
#  of that target. 
.PHONY: all clean debug release


## Default is release build mode
all: release
	
## When in debug mode, don't optimize, and create debug symbols
debug: CXXFLAGS += -O0 -g -D_DEBUG
debug: $(TARGET)
	
## When in release mode, optimize
release: CXXFLAGS += -O3
release: $(TARGET)

## Remove built object files and the main executable
clean:
	$(info ... deleting built object files and executable  ...)
	-rm *.o $(TARGET)

## The main executable depends on all object files of all source files
$(TARGET): $(OBJECTS) Makefile
	$(info ... linking $@ ...)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

## Every object file depends on its source. It may also depend on
#  potentially all header files, and of course the makefile itself.
$(OBJDIR)/%.o: %.c Makefile $(HEADERS)
	$(info ... compiling $@ ...)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $< -o $@

$(OBJDIR)/%.o: %.cc Makefile $(HEADERS)
	$(info ... compiling $@ ...)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $< -o $@

$(OBJDIR)/%.o: %.cpp Makefile $(HEADERS)
	$(info ... compiling $@ ...)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $< -o $@

$(OBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)


