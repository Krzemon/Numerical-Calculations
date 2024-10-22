set term pngcairo
set output '../plots/absorpcja.png'

set xlabel 'X'
set ylabel 'Y'
set title 'Mapa absorpcji'
set xrange [0:200]
set yrange [0:200]

plot '../data/absorpcja.dat' matrix using 2:1:3 with image notitle
