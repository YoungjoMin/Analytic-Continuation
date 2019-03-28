# Makefile for user_programs

all:
        g++ -std=c++11 -g -O2 -I/home/hyunwoo/iRRAM/installed/include -Xlinker -rpath -Xlinker /home/hyunwoo/iRRAM/installed/lib  Compute.cc  -L/home/hyunwoo/iRRAM/installed/lib -liRRAM -lmpfr -lgmp -lm -lpthread -o Compute

