#Adaptation of uC Makefile for native C
#General compiler options and such
CC=gcc
CFLAGS=-O$(OPT) -g -Wall -Wstrict-prototypes
OPT=0#optimization level: s (size), 0 (off), 1, 2, 3

#Source files and output
SRC=main.c#C source files
TARGET=HW1#Output filename
INC=-I./#Directories to search for headers

#Rules

all: run

run: binary
	./$(TARGET).o

binary: $(SRC)
	$(CC) $(CFLAGS) $(INC) $^ -o $(TARGET).o -lm
