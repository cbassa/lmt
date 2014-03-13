# Sample pipeline for converting baseband timeseries

# Make fifos
# The command below will make user specific fifos
mkdir -p /tmp/$USER
if [ ! -e /tmp/$USER/fifo_reader1 ]; then
    mkfifo /tmp/$USER/fifo_reader1
fi
if [ ! -e /tmp/$USER/fifo_channelizer1 ]; then
    mkfifo /tmp/$USER/fifo_channelizer1
fi
if [ ! -e /tmp/$USER/fifo_reader2 ]; then
    mkfifo /tmp/$USER/fifo_reader2
fi
if [ ! -e /tmp/$USER/fifo_channelizer2 ]; then
    mkfifo /tmp/$USER/fifo_channelizer2
fi

# Effelsberg
./dada_reader_nodelay -i /home/bassa/code/c/research/leap/leaptools/data/3C454_EB.dada -o /tmp/$USER/fifo_reader1 -b 64000 -t 7.0 &
./channelizer -i /tmp/$USER/fifo_reader1 -n 160 -o /tmp/$USER/fifo_channelizer1 &

# Westerbork
./dada_reader_nodelay -i /home/bassa/code/c/research/leap/leaptools/data/3C454_WB.dada -o /tmp/$USER/fifo_reader2 -b 64000 -t 0.000486 &
./channelizer -i /tmp/$USER/fifo_reader2 -n 200 -o /tmp/$USER/fifo_channelizer2 &

cat /tmp/$USER/fifo_channelizer1 >out1.dat
cat /tmp/$USER/fifo_channelizer2 >out2.dat