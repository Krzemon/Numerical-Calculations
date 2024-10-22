set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'plots/sigma.png'
set title "Sigma Map"
set xlabel "a index"
set ylabel "c index"
set cblabel "sigma"
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")
set view map
set pm3d at b
splot 'data/sigma.dat' matrix with image
