# Makefile for user_programs

prefix=/usr/local
exec_prefix=/usr/local

CC = gcc -std=c11
CFLAGS = -g -O2 
CPP = gcc -E
CPPFLAGS =  -I${prefix}/include
CXX = g++ -std=c++11
CXXCPP = g++ -E -std=c++11
CXXFLAGS = -g -O2
LDFLAGS = -Xlinker -rpath -Xlinker ${exec_prefix}/lib
LDLIBS =  -L${exec_prefix}/lib -liRRAM -lmpfr -lgmp -lm -lpthread

SOURCE = POWERSERIES Compute 
OBJECT = POWERSERIES.o Compute.o

.PHONY: clean
all:	Link

Link: $(OBJECT)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJECT) $(LDLIBS)

POWERSERIES.o: POWERSERIES.cc POWERSERIES.h
	$(CXX) $(CXXFLAGS) -c $(CPPFLAGS) $(LDFLAGS) POWERSERIES.cc $(LDLIBS)

ANALYTIC.o: ANALYTIC.cc ANALYTIC.h
	$(CXX) $(CXXFLAGS) -c $(CPPFLAGS) $(LDFLAGS) ANALYTIC.cc $(LDLIBS)

Compute.o: Compute.cc
	$(CXX) $(CXXFLAGS) -c $(CPPFLAGS) $(LDFLAGS) Compute.cc $(LDLIBS)

clean:
	rm *.o;rm *.out;
