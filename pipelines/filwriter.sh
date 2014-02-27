# Sample pipeline for creating filterbank files

# Make fifos
if [ ! -e /tmp/$USER/fifo_reader ]; then
    mkfifo /tmp/$USER/fifo_reader
fi
if [ ! -e /tmp/$USER/fifo_channelizer ]; then
    mkfifo /tmp/$USER/fifo_channelizer
fi
if [ ! -e /tmp/$USER/fifo_integrator ]; then
    mkfifo /tmp/$USER/fifo_integrator
fi

# Start reader
# input: input file, block size
dada_reader_nodelay -i input.dada -o /tmp/$USER/fifo_reader -b 64000 &

# Start channelizer
# input: number of channels
channelizer -i /tmp/$USER/fifo_reader -o /tmp/$USER/fifo_channelizer -n 64 &

# Start integrator
# input: number of spectra to average
simple_integrator -i /tmp/$USER/fifo_channelizer -o /tmp/$USER/fifo_integrator -n 100 &

# Start filwriter
# input: output file
filwriter -i /tmp/$USER/fifo_integrator -o out.fil
