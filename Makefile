#       
# Makefile for building a program that calls a routine from
# the CLAPACK Library
#
# En MaC es necesario -I/opt/local/include -L/opt/local/lib/

CC = gcc
CFLAGS =-O3 -Wall 
LLAPK= -lfftw3  -lm

all: cont_complex.exe
	echo $<

cont_complex.exe: cont_complex.o symmetric_curves.o
	$(CC) $(CFLAGS) $^ -o $@  -lfftw3 -lm

cont_complex.o: cont_complex.c parameter.h
	$(CC) $(CFLAGS) $<   -c

symmetric_curves.o: symmetric_curves.c symmetric_curves.h
	$(CC) $(CFLAGS)  $< -c
contb.exe: contb.o
	$(CC) $(CFLAGS) contb.o  -o contb.exe $(LLAPK)

contb.o: contb.c
	                $(CC) $(CFLAGS) contb.c  -c

