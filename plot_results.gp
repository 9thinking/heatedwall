#!/usr/bin/gnuplot
#
# Gnuplot script for visualizing droplet impact simulation results
#

# Set output format
set terminal pngcairo enhanced font 'Arial,12' size 1200,800

# Plot 1: Spreading factor vs time
set output 'spreading_factor.png'
set title 'Droplet Spreading Factor Evolution'
set xlabel 'Non-dimensional Time (t × U₀/D₀)'
set ylabel 'Spreading Factor (D/D₀)'
set grid
set key top right
plot 'quantitative_output.dat' u 1:2 w l lw 2 lc rgb 'blue' title 'D/D₀'

# Plot 2: Kinetic energy vs time
set output 'kinetic_energy.png'
set title 'Droplet Kinetic Energy Evolution'
set xlabel 'Non-dimensional Time (t × U₀/D₀)'
set ylabel 'Kinetic Energy'
set grid
set key top right
plot 'quantitative_output.dat' u 1:6 w l lw 2 lc rgb 'red' title 'KE'

# Plot 3: Center of mass trajectory
set output 'center_of_mass.png'
set title 'Droplet Center of Mass Position'
set xlabel 'Radial Position (r/D₀)'
set ylabel 'Axial Position (z/D₀)'
set grid
set key top right
set size ratio -1
plot 'quantitative_output.dat' u 3:4 w l lw 2 lc rgb 'green' title 'CoM Trajectory'

# Plot 4: Multi-panel overview
set output 'overview.png'
set multiplot layout 2,2 title 'Droplet Impact Simulation Overview'

# Panel 1: Spreading factor
set title 'Spreading Factor'
set xlabel 'Time'
set ylabel 'D/D₀'
set grid
plot 'quantitative_output.dat' u 1:2 w l lw 2 lc rgb 'blue' notitle

# Panel 2: Kinetic energy
set title 'Kinetic Energy'
set xlabel 'Time'
set ylabel 'KE'
set grid
plot 'quantitative_output.dat' u 1:6 w l lw 2 lc rgb 'red' notitle

# Panel 3: Volume conservation
set title 'Volume Conservation'
set xlabel 'Time'
set ylabel 'Volume (V/V₀)'
set grid
plot 'quantitative_output.dat' u 1:($5*6.0) w l lw 2 lc rgb 'orange' notitle

# Panel 4: Center of mass height
set title 'Center of Mass Height'
set xlabel 'Time'
set ylabel 'y/D₀'
set grid
plot 'quantitative_output.dat' u 1:4 w l lw 2 lc rgb 'purple' notitle

unset multiplot

# Reset terminal
set output
print "Plots generated successfully!"
print "  - spreading_factor.png"
print "  - kinetic_energy.png"
print "  - center_of_mass.png"
print "  - overview.png"
