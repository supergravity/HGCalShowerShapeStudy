CXXFLAGS=-g -m64 -O2 -Wall $(shell root-config --cflags ) -I ../utils -I ../BaseParticle -I../Electrons -I ../Muons
LIBS=-g -m64  $(shell root-config --libs)    -lGenVector

LOCAL = .
COMMON = $(LOCAL)
CC:=g++

.PHONY: all clean clear

all: makePlots

VPATH = ../utils ../Muons ../Electrons ../BaseParticle  

makePlots :	main.o makePlots.o 
		$(CC) $^ -o $@ $(LIBS)

main.o :	main.cxx  makePlots.h
		$(CC) -c $(CXXFLAGS) $< -o $@

makePlots.o :    makePlots.cxx   makePlots.h waferGeom.h  CompartmentObject.h
		$(CC) -c $(CXXFLAGS) $< -o $@

EnergyRescalerUpgrade.o :    EnergyRescalerUpgrade.cxx   EnergyRescalerUpgrade.h
		$(CC) -c $(CXXFLAGS) $< -o $@


clean:clear	
clear:	
	rm -f ../*/*.o 
	rm -f makePlots
