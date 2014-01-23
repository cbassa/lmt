# Compiling flags
CFLAGS = -O3

# Linking flags
LFLAGS = -lm -lfftw3f

# Compiler
CC = gcc

reader: reader.o dada.o
	$(CC) -o reader reader.o dada.o $(LFLAGS)

integrator: integrator.o
	$(CC) -o integrator integrator.o $(LFLAGS)

channelizer: channelizer.o
	$(CC) -o channelizer channelizer.o $(LFLAGS)

clean:
	rm -f *.o
	rm -f *~
