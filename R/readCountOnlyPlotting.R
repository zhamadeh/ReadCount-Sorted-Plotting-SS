#Script for plotting read count only and ordering plots by coverage

#User input: files2plot, variable to order libraries by, number of libraries to show
#variable to order libraries by can only be: background.estimate, med.reads.per.MB or perc.coverage
# Run: Rscript readCountOnlyPlotting.R BPR2.0/Input/RData_good/ perc.coverage 10

#args = commandArgs(trailingOnly=TRUE)
library(tidyverse)

plottingReadCounts <- function(files2plot="BPR_output/clean/",feature="perc.coverage",numLibsToShow=10,printFile=paste0("readCounts_",feature,".pdf"),halfHalf=FALSE){
	if (is(files2plot, class.breakpoint)) {
		numplots <- 1
	} else if (is.character(files2plot)) {
		numplots <- length(list.files(files2plot))
	} else {
		stop("Unsupported object class submitted!!!")
	}

	plots <- list()

	df=data.frame()
	for (i in seq_len(numplots)) {
		if (is(files2plot, 'character')) {
			data <- loadFromFiles(list.files(files2plot,full.names = T)[i], check.class=class.breakpoint)[[1]]
		} else if (is(files2plot, class.breakpoint)) {
			data <- files2plot
		} else {
			stop("Only 'BreakPoint' class object can be plotted")
		}
		tmp <- data.frame(data$ID,data$lib.metrics[[feature]])
		df <- rbind(tmp,df)
	}

	if (feature=="background.estimate"){
		df <- df[order(df[,2] ),]
	} else { df <- df[order(-df[,2] ),] }


	for (file in df$data.ID){
		file=paste0(files2plot,file,".RData")
		if (is(files2plot, 'character')) {
			data <- loadFromFiles(file, check.class=class.breakpoint)[[1]]
		} else if (is(files2plot, class.breakpoint)) {
			data <- files2plot
		} else {
			stop("Only 'BreakPoint' class object can be plotted")
		}

		filename <- data$ID
		ptm <- startTimedMessage("Plotting ", filename, " ...")

		bamfile <- data$ID
		reads <- data$fragments
		chroms2plot <- GenomeInfoDb::seqlevels(reads)
		breaks <- data$breaks
		counts <- data$counts
		lib.metrics <- data$lib.metrics
		lib.metrics <- round(lib.metrics, digits = 5)
		lib.metrics <- paste(names(lib.metrics), lib.metrics, sep = '=')
		lib.metrics <- paste(lib.metrics, collapse = "  |  ")

		#Skip chromosomes shorter then 5-times of the bin size 200kb
		if (any(seqlengths(reads) < 200000*5)) {
			message(" Skipping short chromosomes/contigs!")
			keep.chroms <- names(seqlengths(reads)[seqlengths(reads) >= 200000*5])
			reads <- GenomeInfoDb::keepSeqlevels(reads, keep.chroms, pruning.mode = 'coarse')
		}

		binned.data <- unlist(GenomicRanges::tileGenome(seqlengths(reads), tilewidth = 200000))

		#counts overlaps between bins and our reads
		Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='-'])
		Crickreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='+'])
		bothreads <- Watsonreads + Crickreads

		mcols(binned.data)$bothreads <- bothreads
		mcols(binned.data)$Watsonreads <- Watsonreads
		mcols(binned.data)$Crickreads <- Crickreads

		#transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
		cum.seqlengths <- cumsum(as.numeric(seqlengths(binned.data)))
		cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
		#get positions of ends of each chromosome to plot lones between the chromosomes
		if (length(cum.seqlengths) > 1) {
			chr.lines <- data.frame( y=cum.seqlengths[-length(cum.seqlengths)] )
		} else {
			chr.lines <- data.frame( y=0 )
		}
		#get positions of each chromosomes names
		chr.label.pos <- round( cum.seqlengths.0 + (0.5 * seqlengths(binned.data) ) )
		names(chr.label.pos) <- gsub("chr", "", names(chr.label.pos)) #line to add to exclude chr

		#transform chromosome based coordinates into genomewide coordinates
		trans.reads <- transCoord(binned.data)
		trans.breaks <- transCoord(breaks)
		trans.counts <- transCoord(counts)

		dfplot.reads <- as.data.frame(trans.reads)
		dfplot.breaks <- as.data.frame(trans.breaks)
		dfplot.counts <- as.data.frame(trans.counts)

		my_theme <- theme(
			legend.position="none",
			panel.background=element_blank(),
			panel.border=element_blank(),
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			plot.background=element_blank()
		)


		### PLOT READS

		#get midpoint values for each genomic bin
		dfplot.reads$midpoint <- dfplot.reads$start.genome + ( (dfplot.reads$end.genome - dfplot.reads$start.genome) %/% 2 )

		#filter bins with extremely high amount of reads
		Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.999)
		Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.999)
		#set outlier bins to the limit
		dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
		dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier

		#construct ggplot
		dfplot.reads$mCrickreads <- -dfplot.reads$Crickreads
		ggplt1 <-
			ggplot(dfplot.reads) +
			geom_linerange(aes_string(ymin=0, ymax='mCrickreads', x='midpoint'), color="paleturquoise4", size=0.2) +
			geom_linerange(aes_string(ymin=0, ymax='Watsonreads', x='midpoint'), color="sandybrown", size=0.2 ) +
			geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black') + xlab(NULL) +
			ylab("Read counts") +
			xlab("Chromosomes") +
			scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) +
			#theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			my_theme +
			ggtitle(bamfile, subtitle = lib.metrics)

		#if (length(breaks)) {
		#	p <- suppressWarnings( cowplot::plot_grid(ggplt1, ggplt2, ggplt3, ncol=1, align="v", rel_heights = c(3,3,2)) )
		#} else {
		#	p <- suppressWarnings( cowplot::plot_grid(ggplt1, ncol=1, align="v", rel_heights = 3) )
		#}
		p=ggplt1
		plots[[length(plots)+1]] <- p
		stopTimedMessage(ptm)
	}
	if (halfHalf==T){
		top=plots[1:(as.numeric(numLibsToShow)/2)]
		bottom=plots[-1:-(length(plots)-(as.numeric(numLibsToShow))/2)]
		plots = c(top, bottom)
	}
	else { plots <- plots[1:numLibsToShow] }

	if (!is.null(printFile)) {
		message("Printing to PDF ",printFile)

		grDevices::pdf(printFile, width=max(10, length(chroms2plot)), height=2)
		bquiet = lapply(plots, print)
		d <- grDevices::dev.off()
	}
}

transCoord <- function(gr) {
	cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
	cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
	names(cum.seqlengths.0) <- GenomeInfoDb::seqlevels(gr)
	gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
	gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
	return(gr)
}

loadFromFiles <- function(files, check.class=c('GRanges', 'BreakPoint')) {

	# ptm <- startTimedMessage("Loading data from files ...")
	if (is.null(files)) {
		# stopTimedMessage(ptm)
		return(files)
	}
	if (any(! check.class %in% c('GRanges', class.breakpoint))) {
		stop("Argument 'check.class' must contain any combination of c('", paste0(c('GRanges', class.breakpoint), collapse="', '"), "').")
	}
	modellist <- list()
	if (is.character(files)) {
		for (file in files) {
			temp.env <- new.env()
			model <- get(load(file, envir=temp.env), envir=temp.env)
			if (! class(model) %in% check.class) {
				stop("File '", file, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
			}
			modellist[[file]] <- model
		}
	} else if (class(files) %in% check.class) {
		modellist[[1]] <- files
	} else if (is.list(files)) {
		for (file in files) {
			model <- file
			if (! class(model) %in% check.class) {
				stop("List entry '", length(modellist)+1, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
			}
			modellist[[length(modellist)+1]] <- model
		}
		names(modellist) <- names(files)
	} else if (! class(files) %in% check.class) {
		stop("Input does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
	}
	# stopTimedMessage(ptm)
	return(modellist)
}

class.breakpoint <- "BreakPoint"

startTimedMessage <- function(...) {

	x <- paste0(..., collapse='')
	message(x, appendLF=FALSE)
	ptm <- proc.time()
	return(ptm)

}

stopTimedMessage <- function(ptm) {

	time <- proc.time() - ptm
	message(" ", round(time[3],2), "s")

}


### Running
#plottingReadCounts(files2plot=args[1],feature=args[2],numLibsToShow=args[3],halfHalf=args[4])


