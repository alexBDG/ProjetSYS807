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
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:2 w l lt rgb "blue", "Methode_Roe_Correctif_1_M_100_Test_1.dat" u 1:2 w l lt rgb "gray", "Methode_Roe_Correctif_2_M_100_Test_1.dat" u 1:2 w l lt rgb "red", "Methode_Roe_Correctif_3_M_100_Test_1.dat" u 1:2 w l lt rgb "green"#, "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:2 w l lt rgb "black"

# --- pressure
set size 0.4,0.5
set origin 0.0,0.0
unset key
set ytics 0.5
set yrange [0:1.1]
set ylabel 'pressure'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:4 w l lt rgb "blue", "Methode_Roe_Correctif_1_M_100_Test_1.dat" u 1:4 w l lt rgb "gray", "Methode_Roe_Correctif_2_M_100_Test_1.dat" u 1:4 w l lt rgb "red", "Methode_Roe_Correctif_3_M_100_Test_1.dat" u 1:4 w l lt rgb "green"#, "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:4 w l lt rgb "black"

# --- velocity
set size 0.4,0.5
set origin 0.4,0.5
unset key
set ytics 0.75
set yrange [-0.1:1.5]
set ylabel 'velocity'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:3 w l lt rgb "blue", "Methode_Roe_Correctif_1_M_100_Test_1.dat" u 1:3 w l lt rgb "gray", "Methode_Roe_Correctif_2_M_100_Test_1.dat" u 1:3 w l lt rgb "red", "Methode_Roe_Correctif_3_M_100_Test_1.dat" u 1:3 w l lt rgb "green"#, "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:3 w l lt rgb "black"

# --- internal energy
set size 0.4,0.5
set origin 0.4,0.0
unset key
set ytics 1.8
set yrange [1.8:3.6]
set ylabel 'internal energy'
plot "Methode_Roe_Correctif_0_M_100_Test_1.dat" u 1:5 w l lt rgb "blue", "Methode_Roe_Correctif_1_M_100_Test_1.dat" u 1:5 w l lt rgb "gray", "Methode_Roe_Correctif_2_M_100_Test_1.dat" u 1:5 w l lt rgb "red", "Methode_Roe_Correctif_3_M_100_Test_1.dat" u 1:5 w l lt rgb "green"#, "Methode_LaxFriedrichs_M_10000_Test_1.dat" u 1:5 w l lt rgb "black"

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
plot 2 w l lt rgb "blue" title "Sans Correctif", 2 w l lt rgb "gray" title "Correctif HH2", 2 w l lt rgb "red" title "Correctif HH1", 2 w l lt rgb "green" title "Correctif HC"#, 2 lt rgb "white" title " ", 2 w l lt rgb "black" title "Lax-Friedrichs", 2 lt rgb "white" title "M=10000"


unset multiplot

pause -1
