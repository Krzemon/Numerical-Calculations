set terminal pngcairo
set output "./rv_10.png"

set xlabel "X"
set ylabel "Y"

set xrange [-0.01:0.5]
set yrange [-0.01:0.5]
set key off

plot "./rv_10.dat" using 1:2
