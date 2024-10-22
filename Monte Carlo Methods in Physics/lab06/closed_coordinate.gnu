set terminal pngcairo size 800,800
set output 'plot/position_1000.png'

set key on
set key inside top right

set pointsize 1.0

set xrange [-6:6]
set yrange [-6:6]

# Narysuj okrąg o promieniu 5 i środku w (0,0)
set object 1 circle at 0,0 size 5 fillstyle empty border lc rgb "red"

# Narysuj dwa małe okręgi
set object 2 circle at -4.5,0 size 0.5 fillstyle empty border lc rgb "blue"
set object 3 circle at 3,0 size 0.5 fillstyle empty border lc rgb "green"

plot 'data/position_1000.dat' with points lc rgb "black" title 't = 1.0'
