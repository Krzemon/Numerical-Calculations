set term png
set output "../plots/potential_mc.png"
set title "Mapa kolorów potencjału z metody Monte Carlo"
set xlabel "X"
set ylabel "Y"
set palette defined (0 "blue", 1 "red")
set xrange [0:30]
set yrange [0:30]
plot "../data/potential_mc.txt" matrix with image notitle