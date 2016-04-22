#       
# Makefile for building a program that calls a routine from
# the CLAPACK Library
#
# En MaC es necesario -I/opt/local/include -L/opt/local/lib/

CC = gcc
CFLAGS =-O3 -Wall 
LLAPK= -lfftw3  -lm


cont_complex.exe: cont_complex.o 
	$(CC) $(CFLAGS) cont_complex.o  -o cont_complex.exe  -lfftw3 -lm

cont_complex.o: cont_complex.c cont_complex.h parameter.h
	$(CC) $(CFLAGS) cont_complex.c  -c

contb.exe: contb.o
	$(CC) $(CFLAGS) contb.o  -o contb.exe $(LLAPK)

contb.o: contb.c
	                $(CC) $(CFLAGS) contb.c  -c

