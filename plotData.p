# Plotting file
####################################################################################

# New plot
reset

# Set the out file


# General theme (modified from Christoph)
set term post eps color enhanced dashed dashlength 3.3 linewidth 1.5 "Times-Roman" 25
unset title
set key spacing 1.5 samplen 3
set xlabel "t / T"
set ylabel "sigma^2"
set tics scale 1.3
set mxtics 10
set mytics 10
unset log y
unset log x
set size 1.25, 1.25
set origin 0.0, 0.0

outputDat = "output_M2MA10.dat"

set out "M2MA10.eps"

p outputDat u ($2/10):6 with lines linestyle 1 lw 2 lt 1 title 'True dispersion', \
  outputDat u ($2/10):7 with lines linestyle 1 lw 2 lt 2 title 'Prolate Recon.', \
  outputDat u ($2/10):8 with lines linestyle 1 lw 2 lt 3 title 'Oblate Recon.',

# New plot
reset

# Set the out file
set out "M2MA10_cor.eps"

# General theme (modified from Christoph)
set term post eps color enhanced dashed dashlength 3.3 linewidth 1.5 "Times-Roman" 25
unset title
set key spacing 1.5 samplen 3
set xlabel "sigma^2_{true}"
set ylabel "sigma^2_{reconstructed}"
set tics scale 1.3
set mxtics 10
set mytics 10
unset log y
unset log x
set size 1.25, 1.25
set origin 0.0, 0.0

p outputDat u 6:7 lw 2 lt 1 title 'True dispersion', \
  outputDat u 6:8 lw 2 lt 2 title 'Prolate Recon.', \
  x with lines lw 2 lt -1
