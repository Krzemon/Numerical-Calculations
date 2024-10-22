set terminal pngcairo
set output "./plot/Ellipse.png"

set xlabel "Polozenie X"
set ylabel "Polozenie Y"

set xrange [-1:1]
set yrange [-1:1]
set key off

plot "./data/Ellipse.dat" using 1:2 linecolor rgb "black" pointtype 6
