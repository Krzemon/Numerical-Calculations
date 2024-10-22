set terminal pngcairo size 1200,800 enhanced font 'Helvetica,16'
set output "zero_energy.png"

set title "Wykresy unormowanej funkcji falowej w poblizu pierwszego miejsca zerowego"
set xlabel "Polozenie [a_b]"
set ylabel "Psi_N"

set xzeroaxis
set grid
set style line 1 lt 1 lw 2 lc rgb 'blue'
set style line 2 lt 1 lw 2 lc rgb 'red'
set style line 3 lt 1 lw 2 lc rgb 'green'

# Ręcznie dostosuj zakres osi y, na przykład od -0.015 do 0.015
set yrange [*:*] # Automatycznie dostosuj zakres osi y
set xrange [0:1900]
dx = (100.0 / 0.05292) / 100
plot "./data/0_energy_0.txt" using ($0*dx):1 with lines linestyle 1 title "Psi_N(E = 0.560 meV)", \
     "./data/0_energy_p.txt" using ($0*dx):1 with lines linestyle 2 title "Psi_N(E = 0.588 meV)", \
     "./data/0_energy_m.txt" using ($0*dx):1 with lines linestyle 3 title "Psi_N(E = 0.532 meV)"