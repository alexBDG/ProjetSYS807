#!/usr/bin/gnuplot

set terminal qt size 800,500


set xtics 0.5
set xrange [0:1.]
set xlabel 'position'

set multiplot

# --- density
set size 0.4,0.5
set origin 0.0,0.5
unset key
set ytics 0.5
set yrange [0:1.1]
set ylabel 'density'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:2 title "Roe", "Methode_LaxFriedrichs_M_100_Test_1.dat" u 1:2 w p lt rgb "gray" title "Lax-Friedrichs", "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:2 w l lt rgb "black" title "Lax-Friedrichs"

# --- pressure
set size 0.4,0.5
set origin 0.0,0.0
unset key
set ytics 0.5
set yrange [0:1.1]
set ylabel 'pressure'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:4 title "Roe", "Methode_LaxFriedrichs_M_100_Test_1.dat" u 1:4 w p lt rgb "gray" title "Lax-Friedrichs", "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:4 w l lt rgb "black" title "Lax-Friedrichs"

# --- velocity
set size 0.4,0.5
set origin 0.4,0.5
unset key
set ytics 0.75
set yrange [-0.1:1.5]
set ylabel 'velocity'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:3 title "Roe", "Methode_LaxFriedrichs_M_100_Test_1.dat" u 1:3 w p lt rgb "gray" title "Lax-Friedrichs", "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:3 w l lt rgb "black" title "Lax-Friedrichs"

# --- internal energy
set size 0.4,0.5
set origin 0.4,0.0
unset key
set ytics 1.8
set yrange [1.8:3.6]
set ylabel 'internal energy'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:5 title "Roe", "Methode_LaxFriedrichs_M_100_Test_1.dat" u 1:5 w p lt rgb "gray" title "Lax-Friedrichs", "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:5 w l lt rgb "black" title "Lax-Friedrichs"

# --- legend
set size 0.2,0.4
set origin 0.8,0.35
set key box bottom center
set key title "Courbes"
set border 0
unset tics
unset xlabel
unset ylabel
set yrange [0:1]
plot 2 w p title "Roe", 2 lt rgb "white" title "M=100", 2 lt rgb "white" title " ", 2 w p pointtype 2 lt rgb "gray" title "Lax-Friedrichs", 2 lt rgb "white" title "M=100", 2 lt rgb "white" title " ", 2 w l lt rgb "black" title "Lax-Friedrichs", 2 lt rgb "white" title "M=10000"


unset multiplot

pause -1
