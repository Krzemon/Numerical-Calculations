set term png
set output "../plots/absorbed_chains.png"
set title "Mapa kolorów ułamka zaabsorbowanych łańcuchów"
set xlabel "X"
set ylabel "Y"
set palette defined (0 "blue", 1 "red")
set xrange [0:30]
set yrange [0:30]
plot "../data/absorbed_chains.txt" matrix with image notitle