set terminal pngcairo enhanced
set output './plot/plot.png'

set size ratio -1

set xlabel "X"
set ylabel "Y"

plot './data/points_A.dat' using 1:2 title 'K_A' with points pt 7 ps 0.5 lc rgb 'blue', \
     './data/points_B.dat' using 1:2 title 'K_B' with points pt 7 ps 0.5 lc rgb 'red'
