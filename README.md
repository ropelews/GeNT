# GeNT
GeNT is a FORTRAN program that was originally written by Hugh Nicholas at the Pittsburgh Supercomputing Center. It is provided on this github page for people that still find it of value. Here is what the code does:  

       This program computes the cross entropy for groups of sequences
       that have been assigned to groups on the basis of biochemical,
       physiological, or other biological property.  The sequence
       assignments are cross-validated, again by the cross entropy
       measure, to check for problems with the alignment or group
       assignment.  Sequences that were initially identified as
       "unclassified" are compared to all of the groups using position
       specific log-odds scores as described by Henikoff and Henikoff.
       Positions in the aligned sequences that are important for
       determining group membership are identified by having a high
       entropy for the entire alignment and a high entropy for one or
       more specific groups.

 

**To Compile:**

gfortran -o gent gent.for

**To run:**

gent

**Plotting results:** 

The file gnuplot.script worked with gnuplot to plot output (plot.data) into a PNG file (plot.png) from the gent program.
 
