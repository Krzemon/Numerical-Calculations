set terminal gif animate delay 0 optimize size 960, 720
set output './velocity_histogram.gif'
set datafile separator ","

datafile = "build/data/velocity.txt"
set xlabel 'Velocity'
set ylabel 'Frequency'

# Set initial x-axis range
xmin = 0
xmax = 3000
num_bins = 40  # Set a fixed number of bins
binwidth = (xmax - xmin) / num_bins

set style fill solid 0.5
set boxwidth binwidth
set xtic binwidth*5  # 5 jednostek na osi x między podziałkami

stats datafile u 4 nooutput  # Use column 4 for the wypadkowa predkosc
N = STATS_blocks

# Calculate the yrange dynamically based on the maximum frequency across all frames
stats datafile u 4 nooutput
ymax_global = STATS_max
ymin = 0
ymax = 50 
set yrange [ymin:ymax]

do for [i=0:N-1:1] {  # Loop over all frames with step 5

    set title sprintf('Frame %d', i)

    # Set dynamic x-axis range based on the current frame
    set xrange [xmin:xmax]

    # Generate frequency data with adjusted binning
    plot datafile index i u (binwidth*floor($4/binwidth) + binwidth/2):(1.0) smooth frequency with boxes lc rgb "blue" notitle

    pause 0.1  # Optional: Pause for a short time to create the animation effect
}
