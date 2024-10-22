set terminal pngcairo
set output "./plot/Box_Muller.png"

set xlabel "Polozenie X"
set ylabel "Polozenie Y"

set xrange [-4:4]
set yrange [-4:4]
set key off

plot "./data/Box_Muller.dat" using 1:2 linecolor rgb "black" pointtype 6
