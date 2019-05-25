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
plot "Methode_LaxFriedrichs_Correctif_0_M_100_Test_1.dat" u 1:2, "Methode_LaxFriedrichs_Correctif_0_M_1000_Test_1.dat" u 1:2 w l lt -1, "Methode_LaxFriedrichs_Correctif_0_M_10000_Test_1.dat" u 1:2 w l

# --- pressure
set size 0.4,0.5
set origin 0.0,0.0
unset key
set ytics 0.5
set yrange [0:1.1]
set ylabel 'pressure'
plot "Methode_LaxFriedrichs_Correctif_0_M_100_Test_1.dat" u 1:4, "Methode_LaxFriedrichs_Correctif_0_M_1000_Test_1.dat" u 1:4 w l lt -1, "Methode_LaxFriedrichs_Correctif_0_M_10000_Test_1.dat" u 1:4 w l

# --- velocity
set size 0.4,0.5
set origin 0.4,0.5
unset key
set ytics 0.75
set yrange [-0.1:1.5]
set ylabel 'velocity'
plot "Methode_LaxFriedrichs_Correctif_0_M_100_Test_1.dat" u 1:3, "Methode_LaxFriedrichs_Correctif_0_M_1000_Test_1.dat" u 1:3 w l lt -1, "Methode_LaxFriedrichs_Correctif_0_M_10000_Test_1.dat" u 1:3 w l

# --- internal energy
set size 0.4,0.5
set origin 0.4,0.0
unset key
set ytics 1.8
set yrange [1.8:3.6]
set ylabel 'internal energy'
plot "Methode_LaxFriedrichs_Correctif_0_M_100_Test_1.dat" u 1:5, "Methode_LaxFriedrichs_Correctif_0_M_1000_Test_1.dat" u 1:5 w l lt -1, "Methode_LaxFriedrichs_Correctif_0_M_10000_Test_1.dat" u 1:5 w l

# --- legend
set size 0.2,0.4
set origin 0.8,0.45
set key box reverse bottom center
set key title "Précision du maillage"
set border 0
unset tics
unset xlabel
unset ylabel
set yrange [0:1]
plot 2 w p title "M = 100    ", 2 lt -1 title "M = 1000  ", 2 w l title "M = 10000"


unset multiplot

pause -1
