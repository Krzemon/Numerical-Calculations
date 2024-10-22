set terminal gif animate delay 0 optimize size 960, 720
set output './particles.gif'
set datafile separator ","

#system "mkdir -p plot"

datafile = "build/data/position.txt"
set xlabel 'X'
set ylabel 'Y'
set xrange [-10:10]
set yrange [-10:10]

stats datafile u 1 nooutput
N=STATS_blocks

do for [i=0:N-1] {

    set title sprintf('Frame %d', i)

    plot datafile index i u 2:3 with points pointtype 7 pointsize 2.0 lc rgb "black" notitle

    pause 0.1  # Optional: Pause for a short time to create the animation effect
}
