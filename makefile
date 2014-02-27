# Compiling flags
CFLAGS = -O3 -I/home/leap/linux/include

# Linking flags
LFLAGS = -lm -L/home/leap/linux/lib -lfftw3f
LPGFLAGS = -lcpgplot -lpgplot -lX11 -lpng


# Compilers
CC = gcc
F77 = gfortran

all: 
	make reader channelizer integrator filwriter donothing_fb dada_reader_nodelay dechannelizer dada_writer digitizer simple_integrator plotter map

plotter: plotter.o
	$(F77) -o plotter plotter.o $(LFLAGS) $(LPGFLAGS)

map: map.o
	$(CC) -o map map.o $(LFLAGS)

donothing_fb: donothing_fb.o
	$(CC) -o donothing_fb donothing_fb.o $(LFLAGS)

filwriter: filwriter.o dada.o
	$(CC) -o filwriter filwriter.o dada.o $(LFLAGS)

reader: reader.o dada.o lib/delays.o
	$(CC) -o ~/bin/reader reader.o dada.o lib/delays.o $(LFLAGS)

dada_reader_nodelay: dada_reader_nodelay.o dada.o
	$(CC) -o dada_reader_nodelay dada_reader_nodelay.o dada.o $(LFLAGS)

integrator: integrator.o
	$(CC) -o integrator integrator.o $(LFLAGS)

simple_integrator: simple_integrator.o
	$(CC) -o simple_integrator simple_integrator.o $(LFLAGS)

digitizer: digitizer.o
	$(CC) -o digitizer digitizer.o $(LFLAGS)

channelizer: channelizer.o dada.o
	$(CC) -o channelizer channelizer.o dada.o $(LFLAGS)

dada_writer: dada_writer.o
	$(CC) -o dada_writer dada_writer.o $(LFLAGS)

dechannelizer: dechannelizer.o dada.o
	$(CC) -o dechannelizer dechannelizer.o dada.o $(LFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm -f lib/*.o	
	rm -f lib/*~
