set terminal pngcairo size 800,800

set xlabel "Time"
set ylabel "N"
set title "Liczba aktywnych cząstek"
set output "plot/active.png"
plot "data/active.dat" using 1:2 with lines title ' '
