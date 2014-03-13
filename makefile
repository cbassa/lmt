# Compiling flags
CFLAGS = -O3 -I/home/leap/linux/include

# Linking flags
LFLAGS = -lm -L/home/leap/linux/lib -lfftw3f
GSLFLAGS = -lgsl -lgslcblas
LPGFLAGS = -lcpgplot -lpgplot -lX11 -lpng


# Compilers
CC = gcc
F77 = gfortran

all: 
	make reader channelizer filwriter donothing_fb dada_reader_nodelay dechannelizer dada_writer digitizer integrator simple_integrator plotter map correlator

plotter: plotter.o
	$(F77) -o plotter plotter.o $(LFLAGS) $(LPGFLAGS)

skrfi: skrfi.o
	$(CC) -o skrfi skrfi.o $(LFLAGS) $(GSLFLAGS)

correlator: correlator.o
	$(CC) -o correlator correlator.o $(LFLAGS)

map: map.o
	$(CC) -o map map.o $(LFLAGS)

donothing_fb: donothing_fb.o
	$(CC) -o donothing_fb donothing_fb.o $(LFLAGS)

filwriter: filwriter.o dada.o
	$(CC) -o filwriter filwriter.o dada.o $(LFLAGS)

reader: reader.o dada.o lib/delays.o
	$(CC) -o reader reader.o dada.o lib/delays.o $(LFLAGS)

dada_reader_nodelay: dada_reader_nodelay.o dada.o
	$(CC) -o dada_reader_nodelay dada_reader_nodelay.o dada.o $(LFLAGS)

# Note: mpolyco routine calls Tempo2 in order to fold data!
integrator: integrator.o predict.c ppolyco.c mpolyco.c cldj.f djcl.f
	$(CC) -c cldj.f
	$(CC) -c djcl.f
	$(CC) -c predict.c
	$(CC) -o integrator integrator.o predict.o cldj.o djcl.o $(LFLAGS)

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
