set terminal pngcairo size 1200,800 enhanced font 'Helvetica,16'
set output "psi.png"

set title "Wykres kolumny Psi_N w zależności od numeru wiersza"
set xlabel "Energia [meV]"
set ylabel "Psi_N"

set xzeroaxis
set grid
set style line 1 lt 1 lw 2 lc rgb 'blue'

# Dodaj linię o wartości 0
#set arrow from graph 0,first 0 to graph 1,first 0 nohead linestyle 1

#plot "./data/psi_file.txt" using 0:101 with lines linestyle 1 title "Psi_N"

dE_atom = 35.0/5000
plot "./data/psi_file.txt" using ($0 * dE_atom):101 with lines linestyle 1 title "Psi_N"





