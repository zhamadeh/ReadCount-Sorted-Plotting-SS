source("R/breakSeekr.R")
source("R/confidenceInterval.R")
source("R/deltaWCalculator.R")
source("R/GenotypeBreaks.R")
source("R/importReads.R")
source("R/runBreakpointr.R")
source("R/readCountOnlyPlotting.R")
source("R/rwConfig.R")
source("R/exportUCSC.R")
source("R/utils.R")
source("R/summarizeBreaks.R")
source("R/plotting.R")
args = commandArgs(trailingOnly=TRUE)
library(Rsamtools)
install.packages("foreach",repos = "http://cran.us.r-project.org")
library(foreach)

#breakpointR::breakpointr(inputfolder="TestBAMFiles/",outputfolder = "BPR_output",pairedEndReads = T,windowsize=175,binMethod="reads",peakTh=0.3875,min.mapq=7.75,trim=6.5,background=0.15)
#args=c("TestBAMFiles/" ,"BPR_output", "feature" ,"perc.coverage", 10, FALSE)


breakpointr <- function(inputfolder, outputfolder,printer="feature",feature="perc.coverage",numLibsToShow=10,printFile=paste0("readCounts_",feature,".pdf"),halfHalf=FALSE, configfile=NULL,pairedEndReads = F,windowsize=175,binMethod="reads",peakTh=0.3875,min.mapq=7.75,trim=6.5,background=0.15, numCPU=1, reuse.existing.files=FALSE, multi.sizes=NULL, pair2frgm=FALSE, chromosomes=NULL, filtAlt=FALSE, genoT='fisher', zlim=3.291,  minReads=10, maskRegions=NULL, callHotSpots=FALSE, conf=0.99) {

runBreakpointrANDexport <- function(file, datapath, browserpath, config) {
    savename <- file.path(datapath, paste0(basename(file),'.RData'))
    ## Find breakpoints
    if (!file.exists(savename)) {
        tC <- tryCatch({
            breakpoints <- runBreakpointr(bamfile=file, ID=basename(file), pairedEndReads=config[['pairedEndReads']], pair2frgm=config[['pair2frgm']], min.mapq=config[['min.mapq']], filtAlt=config[['filtAlt']], chromosomes=config[['chromosomes']], windowsize=config[['windowsize']], binMethod=config[['binMethod']],  multi.sizes=config[['multi.sizes']], trim=config[['trim']], peakTh=config[['peakTh']], zlim=config[['zlim']], background=config[['background']], genoT=config[['genoT']], minReads=config[['minReads']], maskRegions=maskRegions, conf=config[['conf']])
        }, error = function(err) {
            stop(file,'\n',err)
        })
        save(breakpoints, file=savename)
    } else {
        breakpoints <- get(load(savename))
    }
    ## Write BED file
    savename.breakpoints <- file.path(browserpath,paste0(basename(file), '_breakPoints.bed.gz'))
    if (!file.exists(savename.breakpoints)) {
        breakpointr2UCSC(index=basename(file), outputDirectory=browserpath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks, confidenceIntervals=breakpoints$confint)
    }
}

#========================
### General variables ###
#========================
config <- NULL
if (is.character(configfile)) {
    ## Read config file ##
    errstring <- tryCatch({
        config <- readConfig(configfile)
        errstring <- ''
    }, error = function(err) {
        errstring <- paste0("Could not read configuration file ",configfile)
    })
    if (errstring!='') {
        stop(errstring)
    }
}

## Set createCompositeFile to FALSE [experimental]
createCompositeFile=FALSE

## Put options into list and merge with conf ##
params <- list(numCPU=numCPU, reuse.existing.files=reuse.existing.files, windowsize=windowsize,
               binMethod=binMethod, multi.sizes=multi.sizes, pairedEndReads=pairedEndReads,
               pair2frgm=pair2frgm, chromosomes=chromosomes, min.mapq=min.mapq, filtAlt=filtAlt,
               trim=trim, peakTh=peakTh, zlim=zlim, background=background, genoT=genoT,
               minReads=minReads, createCompositeFile=createCompositeFile, maskRegions=maskRegions,
               callHotSpots=callHotSpots, conf=conf)
config <- c(config, params[setdiff(names(params),names(config))])

## Checks user input ##
if (!config[['pairedEndReads']] & config[['pair2frgm']]) {
    message("Option pair2frgm=TRUE can be used only when pairedEndReads=TRUE\nSetting pair2frgm to FALSE!!!")
    config[['pair2frgm']] <- FALSE
}

if (config[['createCompositeFile']] & config[['callHotSpots']]) {
    message("If createCompositeFile=TRUE breakpoint hotspots can't be called\nSetting callHotSpots to FALSE!!!")
    config[['callHotSpots']] <- FALSE
}

## Helpers ##
windowsize <- config[['windowsize']]

## Set up the directory structure ##
datapath <- file.path(outputfolder,'data')
browserpath <- file.path(outputfolder,'browserfiles')
plotspath <- file.path(outputfolder,'plots')
breakspath <- file.path(outputfolder,'breakpoints')

## Delete old directory if desired ##
if (config[['reuse.existing.files']]==FALSE) {
    if (file.exists(outputfolder)) {
        message("Deleting old directory ",outputfolder)
        unlink(outputfolder, recursive=TRUE)
    }
}

## Creating required directories ##
if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
}
if (!file.exists(datapath)) { dir.create(datapath) }
if (!file.exists(browserpath)) { dir.create(browserpath) }
if (!file.exists(plotspath)) { dir.create(plotspath) }
if (!file.exists(breakspath)) { dir.create(breakspath) }

## Write README file RR
savename <- file.path(outputfolder, 'README.txt')
cat("", file=savename)
cat("BreakpointR outputfolder contains the following folders.\n", file=savename, append=TRUE)
cat("========================================================\n", file=savename, append=TRUE)
cat("- breakpoints: UCSC browser compatible bedgraphs compiling all breakpoints across all single-cell libraries.\n", file=savename, append=TRUE)
cat("               List of all localized breakpoints in all single-cell libraries.\n", file=savename, append=TRUE)
cat("               Locations of breakpoint hotspots if 'callHotSpots=TRUE'.\n", file=savename, append=TRUE)
cat("- browserfiles: UCSC browser formated files with exported reads, deltaWs and breakPoints for every single-cell library.\n", file=savename, append=TRUE)
cat("                UCSC browser formated list of masked genomic regions if set option 'maskRegions'.\n", file=savename, append=TRUE)
cat("- data: RData files storing results of BreakpointR analysis for each single-cell library in an object of class 'BreakPoint'.\n", file=savename, append=TRUE)
cat("- plots: Genome-wide plots for selected chromsosome, genome-wide heatmap of strand states as well as chromosome specific read distribution together with localized breakpoints.\n", file=savename, append=TRUE)

## Make a copy of the config file ##
writeConfig(config = config, configfile=file.path(outputfolder, 'breakpointR.config'))

## Load user defined mask regions ##
if (!is.null(maskRegions)) {
    mask.df <- utils::read.table(maskRegions, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    if (is.null(chromosomes)) {
        maskRegions <- GenomicRanges::GRanges(seqnames=mask.df$V1, IRanges(start=mask.df$V2, end=mask.df$V3))
    } else {
        if (any(mask.df$V1 %in% chromosomes)) {
            mask.df <- mask.df[mask.df$V1 %in% chromosomes,] #select only chromosomes submitted for analysis
            maskRegions <- GenomicRanges::GRanges(seqnames=mask.df$V1, IRanges(start=mask.df$V2, end=mask.df$V3))
        } else {
            warning(paste0('Skipping maskRegions. Chromosomes names conflict: ', mask.df$V1[1], ' not equal ',chromosomes[1]))
            maskRegions <- NULL
        }
    }
}

#=====================================
# Find breakpoints and write BED files
#=====================================
files <- list.files(inputfolder, full.names=TRUE, recursive=FALSE, pattern='.bam$')
if (length(files) == 0) {
    stop("No BAM files present in an 'inputfolder'.")
}

### Run breakpointR ###
if (createCompositeFile) {
    config[['numCPU']] <- 1 #always use only one CPU for composite file analysis
    fragments <- createCompositeFile(file.list=files, chromosomes=config[['chromosomes']], pairedEndReads=config[['pairedEndReads']], pair2frgm=config[['pair2frgm']], min.mapq=config[['min.mapq']], filtAlt=config[['filtAlt']], genoT=config[['genoT']], background=config[['background']])
    runBreakpointrANDexport(file = fragments, datapath = datapath, browserpath = browserpath, config = config)
} else if (numCPU > 1) {
    ## Parallelization ##
    message("Using ",config[['numCPU']]," CPUs")
    cl <- parallel::makeCluster(config[['numCPU']])
    doParallel::registerDoParallel(cl)

    message("Finding breakpoints ...", appendLF=FALSE); ptm <- proc.time()
    temp <- foreach (file = files, .packages=c('breakpointR')) %dopar% {
        runBreakpointrANDexport(file = file, datapath = datapath, browserpath = browserpath, config = config)
    }

    parallel::stopCluster(cl)
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
} else {
    config[['numCPU']] <- 1 #if to use only one CPU or CPU argument not defined
    for (file in files) {
        runBreakpointrANDexport(file = file, datapath = datapath, browserpath = browserpath, config = config)
    }
}

## Write masked regions to BED file
if (!is.null(maskRegions)) {
    ranges2UCSC(gr=maskRegions, index='MaskedRegions', outputDirectory=browserpath, colorRGB='0,0,0')
}

## Compile all breaks using disjoin function
if (createCompositeFile==FALSE) {
    files <- list.files(datapath, pattern=".RData$", full.names=TRUE)

    #breaks.all.files <- GenomicRanges::GRangesList()
    #breaksConfInt.all.files <- GenomicRanges::GRangesList()
    breaks.all.files <- list()
    breaksConfInt.all.files <- list()
    summaryBreaks <- list()
    for (file in files) {
        data <- get(load(file))[c('breaks', 'confint')]
        summaryBreaks[[basename(file)]] <- summarizeBreaks(data)
        breakpoints <- data$breaks
        breaks.confint <- data$confint
        if (length(breakpoints)) {
            suppressWarnings( breaks.all.files[[file]] <- breakpoints ) #TODO check if this can be done without warnings
        }
        if (length(breaks.confint)) {
            suppressWarnings( breaksConfInt.all.files[[file]] <- breaks.confint )
        }
    }

    names(breaks.all.files) <- NULL
    breaks <- do.call(c, breaks.all.files)
    ranges.br <- GenomicRanges::disjoin(breaks) # redefine ranges in df
    hits <- GenomicRanges::countOverlaps(ranges.br, breaks) # counts number of breaks overlapping at each range
    mcols(ranges.br)$hits <- hits # appends hits as a metacolumn in ranges

    names(breaksConfInt.all.files) <- NULL
    breaks.CI <- do.call(c, breaksConfInt.all.files)
    ranges.CI <- GenomicRanges::disjoin(breaks.CI) # redefine ranges in df
    hits <- GenomicRanges::countOverlaps(ranges.CI, breaks) # counts number of breaks overlapping at each range
    mcols(ranges.CI)$hits <- hits # appends hits as a metacolumn in ranges

    ## write a bedgraph of data (overlapping breaks)
    breakpointr2UCSC(index='BreakpointSummary', outputDirectory=breakspath, breaksGraph=ranges.br)
    breakpointr2UCSC(index='BreakpointConfIntSummary', outputDirectory=breakspath, breaksGraph=ranges.CI)
} else {
    files <- list.files(datapath, pattern=".RData$", full.names=TRUE)

    summaryBreaks <- list()
    for (file in files) {
        data <- get(load(file))[c('breaks', 'confint')]
        summaryBreaks[[basename(file)]] <- summarizeBreaks(breakpoints=data)
    }
}

## Write all breakpoints to a file
summaryBreaks.df <- do.call(rbind, summaryBreaks)
summaryBreaks.df$filenames <- rownames(summaryBreaks.df)
write.table(summaryBreaks.df, file = file.path(breakspath, 'breakPointSummary.txt'), quote = FALSE, row.names = FALSE)

## Synchronize read directionality [EXPERIMENTAL]
#files2sync <- list.files(datapath, pattern = ".RData", full.names = TRUE)
#syncReads <- synchronizeReadDir(files2sync)
#breakpointr2UCSC(index="syncReads", outputDirectory=browserpath, fragments=syncReads)

## Search for SCE hotspots ##
if (callHotSpots) {
    hotspots <- hotspotter(gr.list=breaks.all.files, bw = 1000000, pval = 1e-10)
    if (length(hotspots)) {
        ranges2UCSC(gr=hotspots, index='HotSpots', outputDirectory=breakspath, colorRGB='0,0,255')
    }
}

## Plotting
datapath <- file.path(outputfolder,'data')
files2plot <- list.files(datapath, pattern = ".RData", full.names = TRUE)
plotspath <- file.path(outputfolder,'plots')
plotBreakpoints(files2plot=files2plot, file=file.path(plotspath, 'breaksPlot.pdf')) -> beQuiet

if (printer=="feature"){
    files2plot <- file.path(outputfolder,'data/')
    plotspath <- file.path(outputfolder,'plots/')
    plottingReadCounts(files2plot=files2plot,feature=feature,numLibsToShow=numLibsToShow,printFile=paste0(plotspath,printFile),halfHalf=halfHalf)
    }
}

breakpointr(inputfolder=args[1],outputfolder = args[2], printer=args[3],feature=args[4],numLibsToShow = args[5],halfHalf=args[6], pairedEndReads = T,numCPU = 1,windowsize=175,binMethod="reads",peakTh=0.3875,min.mapq=7.75,trim=6.5,background=0.15,multi.sizes=NULL,genoT = "fisher")

#inputfolder=args[1]
#outputfolder = args[2]
#printer=args[3]
#feature=args[4]
#numLibsToShow = args[5]
#halfHalf=args[6]
#pairedEndReads = F
#numCPU = 5
#windowsize=175
#binMethod="reads"
#peakTh=0.3875
#min.mapq=7.75
#trim=6.5
#background=0.15
#multi.sizes=NULL
#genoT = "fisher"

#reuse.existing.files=FALSE
#multi.sizes=NULL
#pair2frgm=FALSE
#chromosomes=NULL
#filtAlt=FALSE
#genoT='fisher'
#zlim=3.291
#minReads=10
#maskRegions=NULL
#callHotSpots=FALSE
#conf=0.99


