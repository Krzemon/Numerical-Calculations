system "mkdir -p plot"

set terminal gif animate delay 0 optimize size 960, 720

set output "plot/grid.gif"

#set view map
set cbrange [0:10]


datafile = "build/data/grid.dat"

stats datafile u 1 nooutput

N=STATS_blocks

do for [i=0:N-2] {
    splot datafile index i u 1:2:4 w pm3d notitle
}
```