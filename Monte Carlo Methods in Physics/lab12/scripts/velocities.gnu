############## Format pliku danych ##############
# x  dens_av  p_av  temp_av  v_av  jx_av  nR_av #
#################################################

# Ustawienia terminala i wyjściowego pliku PNG
set term pngcairo size 800,600 enhanced font 'Verdana,12'

# Wykres rozkładu prędkości w chwili startowej i końcowej
set output 'velocity_distribution_combined.png'
set xlabel 'Pozycja x'
set ylabel 'Prędkość [m/s]'
#set title 'Rozkład prędkości'
set key inside

plot 'nptv_0.dat' using 1:5 title 'Prędkość początkowa' with lines linecolor rgb 'black', \
     'nptv_1k.dat' using 1:5 title 'Prędkość po 1000 iteracji' with lines linecolor rgb 'red', \
     'nptv_5k.dat' using 1:5 title 'Prędkość po 5000 iteracji' with lines linecolor rgb 'blue'
