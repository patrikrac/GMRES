set xtics autofreq
set ytics autofreq
set logscale y
set ticslevel 0
set format y "10^{%T}"
set style data linespoints
set grid
set terminal pdf enhanced
set output "plot_1.pdf"
set title "Convergence 1"
plot "build/convergence_1_20.txt" title "restart=20", "build/convergence_1_40.txt" title "restart=40", "build/convergence_1_80.txt" title "restart=80"
set terminal pdf enhanced
set output "plot_2.pdf"
set title "Convergence 2"
plot "build/convergence_2_20.txt"  title "restart=20", "build/convergence_2_40.txt"  title "restart=40", "build/convergence_2_80.txt"  title "restart=80"
set terminal pdf enhanced
set output "plot_3.pdf"
set title "Convergence 3"
plot "build/convergence_3_20.txt"  title "restart=20", "build/convergence_3_40.txt"  title "restart=40", "build/convergence_3_80.txt"  title "restart=80"
set terminal pdf enhanced
set output "plot_4.pdf"
set title "Convergence 4"
plot "build/convergence_4_20.txt"  title "restart=20", "build/convergence_4_40.txt"  title "restart=40", "build/convergence_4_80.txt"  title "restart=80"
set terminal pdf enhanced
set output "plot_5.pdf"
set title "Convergence 5"
plot "build/convergence_5_20.txt"  title "restart=20", "build/convergence_5_40.txt"  title "restart=40", "build/convergence_5_80.txt"  title "restart=80"