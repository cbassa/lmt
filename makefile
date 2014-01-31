# Compiling flags
CFLAGS = -O3 -I/home/leap/linux/include

# Linking flags
LFLAGS = -lm -L/home/leap/linux/lib -lfftw3f

# Compiler
CC = gcc

all: 
	make reader channelizer integrator filwriter donothing_fb dada_reader_nodelay

donothing_fb: donothing_fb.o
	$(CC) -o donothing_fb donothing_fb.o $(LFLAGS)

filwriter: filwriter.o dada.o
	$(CC) -o filwriter filwriter.o dada.o $(LFLAGS)

reader: reader.o dada.o
	$(CC) -o reader reader.o dada.o $(LFLAGS)

dada_reader_nodelay: dada_reader_nodelay.o dada.o
	$(CC) -o dada_reader_nodelay dada_reader_nodelay.o dada.o $(LFLAGS)

integrator: integrator.o
	$(CC) -o integrator integrator.o $(LFLAGS)

channelizer: channelizer.o dada.o
	$(CC) -o channelizer channelizer.o dada.o $(LFLAGS)

clean:
	rm -f *.o
	rm -f *~
