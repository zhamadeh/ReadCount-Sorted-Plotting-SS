#!/bin/sh

#You can insert this directly in sspipe master script, just  adjust the path to breakpointr.R

####ARGUMENTS####
#1: input folder with BAM files
#2: output folder for BPR output, this is overwrittern and made on the spot
#3: argument is logical for printing full breaksPlot
#4: argument is logical for printing featurer plot
#5: feature to be used in feature plot can either be "perc.coverage", "background.estimate","med.reads.per.MB"
#6: number of files to show in feature file
#7: is FALSE if only top libraries are shown in feature file, if TRUE, shows half top half bottom files
#8: numCPUs to use
#9: to print out SCE summary plots/tables and do metrics quality filtering and AWC blacklisting
#10: metrics file dir for sce summary


Rscript R/breakpointr.R "../Sequencing_data/ALL_BAM_FILES/" "BPR_output" TRUE TRUE "perc.coverage" 10 FALSE 8 TRUE "Input/Metrics/"
