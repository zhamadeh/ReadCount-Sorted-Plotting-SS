#!/bin/sh

#You can insert this directly in sspipe master script, just  adjust the path to breakpointr.R

####ARGUMENTS####
#first is input folder with BAM files
#second if output folder for BPR output, this is overwrittern and made on the spot
#third argument must either be "feature" if you want read-couny only plotting, or NULL if you dont
#fourth argument can either be "perc.coverage", "background.estimate","med.reads.per.MB"
#fifth argument is number of files to show in feature file
#sixth is FALSE if only top libraries are shown in feature file, if TRUE, shows half top half bottom files


Rscript R/breakpointr.R "../NEW_PAPER_OLD_DATA/pe_alignment_pipeline/mdup/" "BPR_output" "feature" "perc.coverage" 10 FALSE
