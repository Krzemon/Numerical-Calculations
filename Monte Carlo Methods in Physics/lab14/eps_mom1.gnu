set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'plots/eps_mom1.png'
set title "Eps Mom1 Map"
set xlabel "a index"
set ylabel "c index"
set cblabel "eps_mom1"
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")
set view map
set pm3d at b
splot 'data/eps_mom1.dat' matrix with image