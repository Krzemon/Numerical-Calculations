############## Format pliku danych ##############
# x  dens_av  p_av  temp_av  v_av  jx_av  nR_av #
#################################################

# Ustawienia terminala i wyjściowego pliku PNG
set term pngcairo size 800,600 enhanced font 'Verdana,12'

# Wykres rozkładu temperatury w czasie wzdłuż kierunku x-owego
set output 'temperature_distribution_5k.png'
set xlabel 'Pozycja x'
set ylabel 'Temperatura [K]'
#set title 'Rozkład temperatury wzdłuż osi x'
plot 'nptv_5k.dat' using 1:4 title 'Temperatura' with lines

# Wykres rozkładu ciśnienia w czasie wzdłuż kierunku x-owego
set output 'pressure_distribution_5k.png'
set xlabel 'Pozycja x'
set ylabel 'Ciśnienie [Pa]'
#set title 'Rozkład ciśnienia wzdłuż osi x'
plot 'nptv_5k.dat' using 1:3 title 'Ciśnienie' with lines

# Wykres rozkładu prędkości w chwili startowej i końcowej
set output 'velocity_distribution_5k.png'
set xlabel 'Pozycja x'
set ylabel 'Prędkość [m/s]'
#set title 'Rozkład prędkości'
plot 'nptv_5k.dat' using 1:5 title 'Prędkość' with lines
