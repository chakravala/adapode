# adapode, Copyright (C) 2013 Michael Reed
# gnuplot script
#set term png size 800, 600

# Parametric Graph
#set output "LorenzParameteric.png"
set xlabel "x1"
set ylabel "x2"
set zlabel "x3"
set title "Lorenz Problem (Michael Reed, MA 448)"
splot "d" using 2:3:4 with lines
pause -1 "Parametric Graph		Press enter to continue"

# Time Graph
#set output "LorenzTimeGraph.png"
set xlabel "Time"
set ylabel "X"
plot "d" using 1:2 with lines, "d" using 1:3 with lines, \
	"d" using 1:4 with lines
pause -1 "Time Graph			Press enter to continue"

# Time Graph (Logarithmic)
#set output "LorenzTimeGraphLogarithmig.png"
set logscale x
plot "d" using 1:2 with lines, "d" using 1:3 with lines, \
	"d" using 1:4 with lines
pause -1 "Time Graph (Logarithmic)	Press enter to continue"
