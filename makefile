# Compiling flags
CFLAGS = -O3

# Linking flags
LFLAGS = -lm -lfftw3f

# Compiler
CC = gcc

all: 
	make reader channelizer integrator

reader: reader.o dada.o
	$(CC) -o reader reader.o dada.o $(LFLAGS)

integrator: integrator.o
	$(CC) -o integrator integrator.o $(LFLAGS)

channelizer: channelizer.o dada.o
	$(CC) -o channelizer channelizer.o dada.o $(LFLAGS)

clean:
	rm -f *.o
	rm -f *~
