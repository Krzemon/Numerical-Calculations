set term pngcairo
set output '../plots/reflektancja.png'

set xlabel 'X'
set ylabel 'Reflektancja'
set title 'Reflektancja'
set xrange [0:200]

plot '../data/reflektancja.dat' using 1 with lines notitle
