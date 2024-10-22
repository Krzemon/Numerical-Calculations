set terminal pngcairo size 800,800
set output 'plot/coordinate.png'

set key on
set key inside top right

set pointsize 0.5
set xrange [-15:15]
set yrange [-15:15]


plot 'data/coordinate_5-0.dat' with points lc rgb "black" title 't = 5.0', \
     'data/coordinate_1-0.dat' with points lc rgb "red" title 't = 1.0', \
     'data/coordinate_0-1.dat' with points lc rgb "blue" title 't = 0.1'
     
     

#set terminal pngcairo enhanced font "arial,10" size 800,600
#set output "particle_positions_plot.png"

#set xlabel "Pozycja x"
#set ylabel "Pozycja y"
#set title "Pozycje cząstek w układzie zamkniętym"

#do for [i=0:999] {
#    plot "pos_".i.".dat" with points pt 7 lc rgb "red" title "Iteracja ".i
#}