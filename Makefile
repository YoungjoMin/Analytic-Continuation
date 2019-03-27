# Makefile for user_programs

prefix=/home/norbert/real/current
exec_prefix=${prefix}

CC = gcc
CFLAGS = -g -O2 -Wall -Wno-implicit -Wmissing-prototypes
CPP = gcc -E
CPPFLAGS =  -I${prefix}/include
CXX = g++
CXXCPP = g++ -E
CXXFLAGS = -g -O2 -Wall
LDFLAGS = -Xlinker -rpath -Xlinker ${exec_prefix}/lib64
LDLIBS =  -L${exec_prefix}/lib64 -liRRAM -lmpfr -lm -lgmp -lpthread

SOURCE = Poly Compute ANAYTIC
OBJECT = Poly.o Compute.o ANALYTIC.o

.PHONY: clean
all:	Link

Link: $(OBJECT)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJECT) $(LDLIBS)

Poly.o: Poly.cc Poly.h
	$(CXX) $(CXXFLAGS) -c $(CPPFLAGS) $(LDFLAGS) Poly.cc $(LDLIBS)

Compute.o: Compute.cc
	$(CXX) $(CXXFLAGS) -c $(CPPFLAGS) $(LDFLAGS) Compute.cc $(LDLIBS)

clean:
	for i in *.o ; do rm -vf `basename $$i .c` ; done;
