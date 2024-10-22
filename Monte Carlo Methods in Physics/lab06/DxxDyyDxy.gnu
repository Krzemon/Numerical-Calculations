set xlabel "Time"
set ylabel "Diffusion Coefficients"
set title "Diffusion Coefficients over Time"
set terminal pngcairo size 800,800
set output "plot/diffusion_coefficients.png"
plot "data/coefficients.dat" using 1:2 with lines title 'D_{xx}', \
     "data/coefficients.dat" using 1:3 with lines title 'D_{xy}', \
     "data/coefficients.dat" using 1:4 with lines title 'D_{yy}'
