set terminal pngcairo size 1280,720 enhanced font 'Arial, 14'
set output 'hist_02k.png'

# Customize these variables as needed
file_name = "hist_02k.dat"
n_mix = 1 # Number of particle types
nhist = 50 # Number of histogram bins

# set title "Velocity Histogram"
set xlabel "Velocity"
set ylabel "Frequency"
set grid

# Define colors for different particle types
colors = "red blue green orange purple brown pink gray cyan magenta"

set style data histograms
set style fill solid border -1
set boxwidth 0.9 relative

# Create separate plot commands for each type of particle
plot for [ic=1:n_mix] file_name using (column(1)):(column(ic+1)) with boxes lc rgb word(colors, ic) title sprintf("Numeric Histogram"), \
     for [ic=1:n_mix] file_name using (column(1)):(column(n_mix+ic+1)) with lines lw 2 lc rgb word(colors, ic) title sprintf("Theoretical Histogram")