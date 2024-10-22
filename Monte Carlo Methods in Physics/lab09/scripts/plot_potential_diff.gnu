set term png
set output "../plots/potential_diff.png"
set title "Mapa kolorów różnicy potencjałów MC i SOR"
set xlabel "X"
set ylabel "Y"
set palette defined (-0.5 "blue", 0 "white", 0.5 "red")
set xrange [0:30]
set yrange [0:30]
plot "../data/potential_diff.txt" matrix with image notitle