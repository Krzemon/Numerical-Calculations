set term png
set output "../plots/potential_sor.png"
set title "Mapa kolorów potencjału SOR"
set xlabel "X"
set ylabel "Y"
set palette defined (-1 "blue", 0 "white", 1 "red")
set xrange [0:30]
set yrange [0:30]
plot "../data/potential_sor.txt" matrix with image notitle