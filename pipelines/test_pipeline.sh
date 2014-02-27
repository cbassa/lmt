# Make fifos
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

# Start reader
./dada_reader_nodelay -i test.dada -o /tmp/$USER/fifo_reader -b 64000 &

# Start channelizer
./channelizer -i /tmp/$USER/fifo_reader -o /tmp/$USER/fifo_channelizer -n 64 &

# Start dechannelizer
./dechannelizer -i /tmp/$USER/fifo_channelizer -o out.dat 
#/tmp/$USER/fifo_dechannelizer &

# Start writer
#./dada_writer -i /tmp/$USER/fifo_dechannelizer -o out.dat -b 64000

# Start integrator
#./integrator -i /tmp/fifo_channelizer -o /tmp/fifo_integrator -n 50 &

# Start filwriter
#./filwriter -i /tmp/fifo_integrator -o test.fil
