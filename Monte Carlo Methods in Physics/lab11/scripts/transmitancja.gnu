set term pngcairo
set output '../plots/transmitancja.png'

set xlabel 'X'
set ylabel 'Transmitancja'
set title 'Transmitancja'
set xrange [0:200]

plot '../data/transmitancja.dat' using 1 with lines notitle
