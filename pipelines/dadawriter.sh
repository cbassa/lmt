# Sample pipeline for converting baseband timeseries

# Make fifos
# The command below will make user specific fifos
mkdir -p /tmp/$USER
if [ ! -e /tmp/$USER/fifo_reader ]; then
    mkfifo /tmp/$USER/fifo_reader
fi
if [ ! -e /tmp/$USER/fifo_channelizer ]; then
    mkfifo /tmp/$USER/fifo_channelizer
fi
if [ ! -e /tmp/$USER/fifo_dechannelizer ]; then
    mkfifo /tmp/$USER/fifo_dechannelizer
fi
if [ ! -e /tmp/$USER/fifo_digitizer ]; then
    mkfifo /tmp/$USER/fifo_digitizer
fi
if [ ! -e /tmp/$USER/fifo_integrator ]; then
    mkfifo /tmp/$USER/fifo_integrator
fi

# Start reader
# input: input file, block size
dada_reader_nodelay -i input.dada -o /tmp/$USER/fifo_reader -b 64000 &

# Start channelizer
# input: number of channels
channelizer -i /tmp/$USER/fifo_reader -n 80 -o /tmp/$USER/fifo_channelizer &

# Start dechannelizer
# input: add -r flag for real output, remove -r flag for complex output
dechannelizer -i /tmp/$USER/fifo_channelizer -o /tmp/$USER/fifo_dechannelizer -r &

# Start digitizer
# input: -s scale -O offset, default are scale=1, offset=0
digitizer -i /tmp/$USER/fifo_dechannelizer -o /tmp/$USER/fifo_digitizer &

# Start writer
# input: output file
dada_writer -i /tmp/$USER/fifo_digitizer -o out.dada

