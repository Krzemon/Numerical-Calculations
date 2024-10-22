set terminal pngcairo
set output "plot/B.png"

set xlabel "n"
set ylabel "common area"

set format x "10^{%L}"

set logscale x 10

set xrange [1:10000000]
set yrange [3:3.5]

plot "./data/B.dat" using (10**$1):2:3 with yerrorbars linecolor rgb "blue" linewidth 2 pointtype 7 pointsize 1 title ""
