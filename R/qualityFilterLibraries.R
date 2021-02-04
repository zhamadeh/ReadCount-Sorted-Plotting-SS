

qualityFilterLibraries <- function(datapath,metricsfileDir,filteredDatapath){

	metrics <- data.frame()

	for (file in list.files(metricsfileDir,full.names = T)){
		met <- read.table(file,header = T,fill=T,row.names = NULL)
		if ("Library" %in% colnames(met) && "Quality" %in% colnames(met)){
			met <- met %>% select(Library,Quality)
			met$file = basename(file)
			}
		met$Library=paste0(met$Library,".trimmed.mdup.bam")
		metrics <- rbind(met,metrics)
	}

	message("\nThere are ",length(metrics$Library)," libraries in the base metrics file with ",ncol(metrics)," features\n")


	#filter desired quality of library
	goodLibraries <- filter(metrics, metrics$Quality=="good")
	goodLibrariesFound <- intersect(paste0(goodLibraries[,1],".RData"),list.files(datapath))

	#print out desired files that weren't found
	if (length(setdiff(goodLibraries$Library,goodLibrariesFound) != 0)){
		message("The following libraries are missing ...\n")
		if (length(setdiff(goodLibraries$Library,goodLibrariesFound) > 5)){
			message(length(setdiff(goodLibraries$Library,goodLibrariesFound))," ... too many to print ... saving to Output/missingLibraries.txt\n")
			write.table(setdiff(goodLibraries$Library,goodLibrariesFound),"Output/missingLibraries.txt",sep="\t",quote=F,row.names = F,col.names = F)
		} else {print(setdiff(paste0(goodLibraries$Library,".RData"),goodLibrariesFound)) }
	}

	#Move good libraries to output directory
	list.of.good.files <- paste0(datapath,'/',goodLibrariesFound)
	file.copy(list.of.good.files, filteredDatapath)
	write.table(goodLibrariesFound,"Output/goodLibrariesFound.txt",sep="\t",quote=F,row.names = F,col.names = F)

	#summary quality stats
	goodLibrariesFound <- as.data.frame(goodLibrariesFound)
	goodLibrariesFound$library = goodLibrariesFound$goodLibrariesFound

	suppressWarnings(goodLibrariesFound <- goodLibrariesFound %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
	goodLibrariesFound$gene <- "gene"
	for (row in 1:nrow(goodLibrariesFound)){

		for (letter in c("a","b","c","d","e","f")){
			if (is.na(goodLibrariesFound[row,letter])!=T){
				if (goodLibrariesFound[row,letter]=="WT" | goodLibrariesFound[row,letter]=="wt"){
					goodLibrariesFound[row,"gene"]="WT"
				}
				else if (goodLibrariesFound[row,letter]=="blm" | goodLibrariesFound[row,letter]=="BLM" ) {
					goodLibrariesFound[row,"gene"]="BLM"
				}

				else if (goodLibrariesFound[row,letter]=="RECQL5" | goodLibrariesFound[row,letter]=="recql5" | goodLibrariesFound[row,letter]=="RECQ5" | goodLibrariesFound[row,letter]=="recq5" ) {
					if (goodLibrariesFound[row,"gene"]=="BLM"){
						goodLibrariesFound[row,"gene"]="BLM/RECQL5"
					}
					else{
						goodLibrariesFound[row,"gene"]="RECQL5"
					}
				}
			}
		}

	}

	goodLibrariesFound <- select(goodLibrariesFound,-c(a,b,c,d,e,f))
	suppressMessages(summaryStats <- group_by(goodLibrariesFound,gene) %>% dplyr::summarize(n()))
	summaryStats$perc <- paste0(round((summaryStats$`n()`/sum(summaryStats$`n()`))*100,2),"%")

	message("Those ",nrow(goodLibrariesFound)," libraries have been added to ",filteredDatapath,"\n\n")

	print(as.data.frame(summaryStats))
}


frequencyFilterBreakpoints <- function(summaryBreaks.df, blacklist,filteredDatapath,cleanDatapath,filterFrequency){

	summaryBreaks.df<- GRanges(summaryBreaks.df)

	centromeres <- read.table(blacklist,header=F) #%>% select(-c(V4))
	centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
	centroGRange <- GRanges(centromeres)

	suppressWarnings(breakGRanges <-summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])
	breakGRanges$freq = 0


	for (row in 1:length(breakGRanges)){
		tmp = breakGRanges[row,]
		breakGRanges[row,]$freq = countOverlaps(tmp,breakGRanges, maxgap = 1000000)
	}

	breakGRanges$allele_freq=breakGRanges$freq/(length(list.files(filteredDatapath)))


	inversions <- breakGRanges[breakGRanges$allele_freq>filterFrequency,]
	breaks <- breakGRanges[breakGRanges$allele_freq<filterFrequency,]
	breaks <- as.data.frame(breaks)
	breaks <- filter(breaks,width<1000000)

	write.table(breaks,"Output/filteredBreakpoints.txt",row.names = F,col.names = T,quote = F,sep="\t")

	breaks$filenames <- tools::file_path_sans_ext(breaks$filenames)
	breaks$filenames <- as.factor(breaks$filenames)

	counter=0
	for (file in list.files(filteredDatapath,full.names = T)){
		counter=counter+1
		message("Working on ",basename(file), " ... ", (counter/length(list.files(filteredDatapath,full.names = T)))*100,"%")
		tmp <- get(load(file))

		seqinfo <- tmp$breaks@seqinfo
		seqnameLevels <- levels(tmp$breaks@seqnames)


		tmp$breaks <-GRanges(filter(breaks,filenames==paste0(tmp$ID,".RData")) %>% select(c(seqnames,start,end,width,strand,genoT,deltaW)))
		tmp$breaks@seqinfo <- seqinfo
		levels(tmp$breaks@seqnames) <- seqnameLevels

		tmp$confint<-GRanges(filter(breaks,filenames==paste0(tmp$ID,".RData")) %>% select(c(seqnames,CI.start,CI.end,width,strand,genoT,deltaW)))
		tmp$confint@seqinfo <- seqinfo
		levels(tmp$confint@seqnames) <- seqnameLevels

		save(tmp, file=paste0(cleanDatapath,"/",basename(file)))
	}

}
