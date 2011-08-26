# Makefile for anmodel

SHELL = /bin/sh
#EXEPATH = ./
EXEPATH = ./
#### system configuration section
# VPATH = ./src
# vpath %.h ./src

#CC = gcc
CC = g++
#CC=c++
#CFLAGS = -ggdb
#CFLAGS = -I./src
#CFLAGS = -O3
#CFLAGS = -O2
#CFLAGS = -O

OBJdependent = hc.hpp synapse.hpp cmpa.hpp filters.hpp
OBJanmodel = cmpa.o synapse.o filters.o anmodel.o complex.o hc.o
#### various target
#
anmodel : $(OBJanmodel)
	$(CC) $(CFLAGS) -lm -o $(EXEPATH)$@ $(OBJanmodel)
model.o : model.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c model.cpp
runmodel.o : runmodel.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c runmodel.cpp
anmodel.o : anmodel.cpp $(OBJdependennt)
	$(CC) $(CFLAGS) -c anmodel.cpp                                                        
synapse.o : synapse.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c synapse.cpp
cmpa.o 	: cmpa.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c cmpa.cpp
filters.o : filters.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c filters.cpp
hc.o	: hc.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c hc.cpp
complex.o : complex.cpp $(OBJdependent)
	$(CC) $(CFLAGS) -c complex.cpp
