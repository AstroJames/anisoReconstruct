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

outputDat = "output_M2MA0.1.dat"

set out "M2MA01_dispersion.eps"

p outputDat u ($2/10):6 with lines linestyle 1 lw 2 lt 1 title 'True dispersion', \
  outputDat u ($2/10):7 with lines linestyle 1 lw 2 lt 2 title 'Prolate Recon.', \
  outputDat u ($2/10):8 with lines linestyle 1 lw 2 lt 3 title 'Oblate Recon.', \
  outputDat u ($2/10):9 with lines linestyle 1 lw 2 lt 4 title 'Average.'
# New plot
reset

# Set the out file
set out "M2MA01_skewness.eps"

# General theme (modified from Christoph)
set term post eps color enhanced dashed dashlength 3.3 linewidth 1.5 "Times-Roman" 25
unset title
set key spacing 1.5 samplen 3
set xlabel "t/T"
set ylabel "Skewness"
set tics scale 1.3
set mxtics 10
set mytics 10
unset log y
unset log x
set size 1.25, 1.25
set origin 0.0, 0.0

p outputDat u ($2/10):15 with lines linestyle 1 lw 2 lt 1 title 'True skewness', \
  outputDat u ($2/10):16 with lines linestyle 1 lw 2 lt 2 title 'Prolate Recon.', \
  outputDat u ($2/10):17 with lines linestyle 1 lw 2 lt 3 title 'Oblate Recon.', \
  outputDat u ($2/10):18 with lines linestyle 1 lw 2 lt 4 title 'Average.', \
# New plot
reset

# Set the out file
set out "M2MA01_kurt.eps"

# General theme (modified from Christoph)
set term post eps color enhanced dashed dashlength 3.3 linewidth 1.5 "Times-Roman" 25
unset title
set key spacing 1.5 samplen 3
set xlabel "t/T"
set ylabel "Kurtosis"
set tics scale 1.3
set mxtics 10
set mytics 10
unset log y
unset log x
set size 1.25, 1.25
set origin 0.0, 0.0

p outputDat u ($2/10):23 with lines linestyle 1 lw 2 lt 1 title 'True kurtosis', \
  outputDat u ($2/10):24 with lines linestyle 1 lw 2 lt 2 title 'Prolate Recon.', \
  outputDat u ($2/10):25 with lines linestyle 1 lw 2 lt 3 title 'Oblate Recon.', \
  outputDat u ($2/10):26 with lines linestyle 1 lw 2 lt 4 title 'Average.'
