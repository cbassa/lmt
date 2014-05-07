# Sample pipeline for converting baseband timeseries

# Make fifos
# The command below will make user specific fifos
mkdir -p /tmp/$USER
rm -f /tmp/$USER/fifo_reader
mkfifo /tmp/$USER/fifo_reader
rm -f /tmp/$USER/fifo_channelizer
mkfifo /tmp/$USER/fifo_channelizer
rm -f /tmp/$USER/fifo_dechannelizer
mkfifo /tmp/$USER/fifo_dechannelizer
rm -f /tmp/$USER/fifo_digitizer
mkfifo /tmp/$USER/fifo_digitizer

# Start reader
# input: input file, block size
../dada_reader_nodelay -i /home/sanidas/testcorr/3C454_WB.dada -o /tmp/$USER/fifo_reader -b 64000 &

# Start channelizer
# input: number of channels
../channelizer -i /tmp/$USER/fifo_reader -n 100 -o /tmp/$USER/fifo_channelizer &

# Start dechannelizer
# input: add -r flag for real output, remove -r flag for complex output
../dechannelizer -i /tmp/$USER/fifo_channelizer -o /tmp/$USER/fifo_dechannelizer &

# Start digitizer
# input: -s scale -O offset, default are scale=1, offset=0
../digitizer -i /tmp/$USER/fifo_dechannelizer -o /tmp/$USER/fifo_digitizer &

# Start writer
# input: output file
../dada_writer -i /tmp/$USER/fifo_digitizer -o out.dada -s ../Dada_header.txt