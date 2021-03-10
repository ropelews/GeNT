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

**Sample Group File**

Below is a sample group file. Sequences are from uniprot database:

       Problem  pa2_small
       Title  PA2 Sequences.
       Group  One
       sp|P00593|PA21B_BOVIN
       sp|P04417|PA21B_AGKHA
       sp|P14418|PA21B_AGKHP
       sp|P81243|PA21B_BOTJA
       sp|P0CAS3|PA212_CRODR
       sp|P0CAS4|PA213_CRODR
       sp|P0CAS5|PA215_CRODU
       sp|P0CAS6|PA216_CRODU
       sp|P0CAS7|PA217_CRODU
       //
       Group  Two
       sp|P00606|PA20_BUNMU
       sp|Q8UW08|PA20_LAPHA
       sp|P08873|PA20_NOTSC
       sp|P20254|PA20_PSEAU
       sp|Q9PUH4|PA210_AUSSU
       sp|Q9PUH3|PA211_AUSSU
       sp|Q9PUH2|PA212_AUSSU
       sp|Q9PUH1|PA213_AUSSU
       sp|Q9PUH0|PA214_AUSSU
       sp|Q9PUG9|PA215_AUSSU
       sp|Q9PUG8|PA216_AUSSU
       sp|Q9PUG7|PA217_AUSSU
       sp|P81236|PA21B_ACAAN
       sp|P59359|PA21B_AUSSU
       sp|Q8QFW4|PA21B_BUNCE
       sp|Q90WA7|PA21B_BUNFA
       sp|P00617|PA21B_BUNMU
       //




**To Compile:**

gfortran -o gent gent.for

**To run:**

gent

**Plotting results:** 

The file gnuplot.script worked with gnuplot to plot program output (plot.data) into a PNG file (plot.png).
 
