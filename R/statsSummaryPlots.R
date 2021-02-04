
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(tidyverse)))


plottingSummary <- function(plotspath,breakspath,metricsfileDir){
	breakpointSummary <- read.table(file.path(breakspath, 'breakPointSummary.txt'),header=T) %>% select(-c(CI.start,CI.end,genoT))

	for (i in seq(1:5)){
		breakpointSummary$filenames <- tools::file_path_sans_ext(breakpointSummary$filenames)
	}

	breakpointSummary$gene <- "gene"
	breakpointSummary$Library <- breakpointSummary$filenames
	suppressWarnings(breakpointSummary <- breakpointSummary %>% separate(Library, c("a","b","c","d","e","f"), "[_-]+"))

	for (row in 1:nrow(breakpointSummary)){
		for (letter in c("a","b","c","d","e","f")){
			#cat(breakpointSummary[1,letter])

			if (breakpointSummary[row,letter]=="WT" | breakpointSummary[row,letter]=="wt"){
				breakpointSummary[row,"gene"]="WT"
			}
			else if (breakpointSummary[row,letter]=="blm" | breakpointSummary[row,letter]=="BLM" ) {
				breakpointSummary[row,"gene"]="BLM"
			}

			else if (breakpointSummary[row,letter]=="RECQL5" | breakpointSummary[row,letter]=="recql5" | breakpointSummary[row,letter]=="RECQ5" | breakpointSummary[row,letter]=="recq5" ) {
				if (breakpointSummary[row,"gene"]=="BLM"){
					breakpointSummary[row,"gene"]="BLM/RECQL5"
				}
				else{
					breakpointSummary[row,"gene"]="RECQL5"
				}
			}
		}

	}

	breakpointSummary <- select(breakpointSummary,-c(a,b,c,d,e,f))

	files <- list.files(metricsfileDir, full.names=T)
	metrics <- data.frame()


	if (length(files)!=0){
		#file=files[4]

		for (file in files){
			met <- read.table(file,header=T,fill=T)
			if ("Reads_aligned_postfiltering" %in% colnames(met)){
				met <- met %>% select(Library, Reads_aligned_postfiltering,Background,Reads_per_Mb)
				names(met)[names(met) == "Reads_aligned_postfiltering"] <- "Reads"
				met$Reads <- as.numeric(gsub("[^0-9.-]", "", met$Reads))
				met[,'Background'] <- as.character(met[,'Background'])
				met[,'Background'] <- str_remove(met[,'Background'],"%")
				met[,'Background']<-as.numeric(met[,'Background'])
			}
			else if ("Postfiltering_reads_aligned" %in% colnames(met)){
				met <- met %>% select(Library, Postfiltering_reads_aligned , Background,Reads_per_Mb)
				names(met)[names(met) == "Postfiltering_reads_aligned"] <- "Reads"
				met$Reads <- as.numeric(gsub("[^0-9.-]", "", met$Reads))
				met[,'Background'] <- as.character(met[,'Background'])
				met[,'Background'] <- str_remove(met[,'Background'],"%")
				met[,'Background']<-as.numeric(met[,'Background'])
			}
			metrics <- rbind(met,metrics)

		}
	}

	breakpointSummary <- merge(breakpointSummary,metrics,by.x="filenames",by.y='Library')

	breakpointSummary$gene <- as.factor(breakpointSummary$gene)
	breakpointSummary$filenames <- as.factor(breakpointSummary$filenames)

	#sces per chroomosome
	suppressMessages(all <- as.data.frame(breakpointSummary %>% group_by(seqnames) %>% dplyr::summarize(ALL=n())))
	suppressMessages(b <- as.data.frame(breakpointSummary %>% filter(gene=="BLM")  %>% group_by(seqnames) %>% dplyr::summarize(BLM=n())))
	suppressMessages(br <- as.data.frame(breakpointSummary %>% filter(gene=="BLM/RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize("BLM/RECQL5"= n())))
	suppressMessages(r <- as.data.frame(breakpointSummary %>% filter(gene=="RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize(RECQL5=n())))
	suppressMessages(w <- as.data.frame(breakpointSummary %>% filter(gene=="WT")  %>% group_by(seqnames)%>% dplyr::summarize(WT=n())))
	byChr <- merge(b,br,by="seqnames")
	byChr <- merge(byChr,r,by="seqnames",all=T)
	byChr <- merge(byChr,w,by="seqnames",all=T)
	byChr <- merge(byChr,all,by="seqnames",all=T)
	byChr[is.na(byChr)] <- 0
	write.table(byChr,"Output/Tables/perChrom.txt",quote=F,row.names = F,col.names = T,sep="\t")

	breakpointSummary$Library=breakpointSummary$filenames
	numOfLibsPerGene <- data.frame(gene=character(),n=numeric())
	b <- (as.data.frame(breakpointSummary %>% filter(gene=="BLM")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
	numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM",n=as.numeric(length(levels(b$Library))))
	br <- (as.data.frame(breakpointSummary %>% filter(gene=="BLM/RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
	numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM.RECQL5",n=as.numeric(length(levels(br$Library))))
	r <- (as.data.frame(breakpointSummary %>% filter(gene=="RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
	numOfLibsPerGene= add_row(numOfLibsPerGene, gene="RECQL5",n=as.numeric(length(levels(r$Library))))
	w <- (as.data.frame(breakpointSummary %>% filter(gene=="WT")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
	numOfLibsPerGene= add_row(numOfLibsPerGene, gene="WT",n=as.numeric(length(levels(w$Library))))
	numOfLibsPerGene= add_row(numOfLibsPerGene, gene="ALL",n=as.numeric(length(levels(breakpointSummary$Library))))
	write.table(numOfLibsPerGene,"Output/Tables/numOfLibsPerGene.txt",quote=F,row.names = F,col.names = F,sep="\t")


	numOfLibsPerGene=read.table("Output/Tables/numOfLibsPerGene.txt",header=F)
	lengths <- read.table("Output/Tables/chrLengths.txt",header=T)
	byChr<- read.table("Output/Tables/perChrom.txt",header=T)

	for (i in c("BLM","BLM.RECQL5","RECQL5","WT","ALL")){
		#print(byChr[,i])
		for (j in 1:nrow(numOfLibsPerGene)){
			if (numOfLibsPerGene[j,1]==i){
				#print(numOfLibsPerGene[j,2])
				byChr[,i] = byChr[,i]/ numOfLibsPerGene[j,2]
			}
		}
	}

	scePerChrPerGeneVsLength <- merge(lengths,byChr,by.x="Chromosome",by.y="seqnames")
	tidy <- gather(scePerChrPerGeneVsLength, gene,sce_per_chr,BLM:ALL)
	all <- filter(tidy,gene=="ALL")
	tidy <- filter(tidy,gene!="ALL")

	suppressMessages(ggplot(tidy) + geom_point(aes(Length,sce_per_chr,group=gene,color=gene), na.rm=TRUE)+
					 	geom_smooth(aes(Length,sce_per_chr,group=gene,color=gene),se=F,method="lm", na.rm=TRUE)+
					 	theme_classic(base_size = 18) +
					 	ylab("SCEs/chr/lib")+
					 	xlab("Chromosome Length") +
					 	ggsave("Output/Plots/scePerChrPerGeneVsLength.png"))

	suppressMessages(ggplot(tidy) + geom_point(aes(Length,sce_per_chr,group=gene,color=gene), na.rm=TRUE)+
					 	geom_smooth(all,mapping=aes(Length,sce_per_chr),color="black",method="lm",size=2, na.rm=TRUE)+
					 	theme_classic(base_size = 18) +
					 	ylab("SCEs/chr/lib")+
					 	xlab("Chromosome Length") +
					 	ggsave("Output/Plots/scePerChrPerGeneVsLength_ALL.png"))




	suppressMessages(test<-as.data.frame(breakpointSummary %>%
										 	group_by(Library) %>%
										 	dplyr::summarize(n())))
	test$gene <- "gene"
	test$library <- test$Library
	suppressWarnings(test <- test %>% separate(Library, c("a","b","c","d","e","f"), "[_-]+"))

	for (row in 1:nrow(test)){
		for (letter in c("a","b","c","d","e","f")){
			#print(test[1,letter])

			if (test[row,letter]=="WT" | test[row,letter]=="wt"){
				test[row,"gene"]="WT"
			}
			else if (test[row,letter]=="blm" | test[row,letter]=="BLM" ) {
				test[row,"gene"]="BLM"
			}

			else if (test[row,letter]=="RECQL5" | test[row,letter]=="recql5" | test[row,letter]=="RECQ5" | test[row,letter]=="recq5" ) {
				if (test[row,"gene"]=="BLM"){
					test[row,"gene"]="BLM/RECQL5"
				}
				else{
					test[row,"gene"]="RECQL5"
				}
			}
		}

	}

	test <- select(test,c("n()","gene"))
	test<-dplyr::rename(test,c("sces"="n()"))

	ggplot(test) + geom_jitter(aes(gene,sces, color=gene))+ geom_boxplot(aes(gene,sces),width=0.1,coef = 5) +
		theme_classic()+
		theme(text=element_text(size=15)) +
		ggsave("Output/Plots/SCEperGene.png")


	#my_comparisons <- list( c("WT", "RECQL5"), c("WT", "BLM/RECQL5"), c("WT", "BLM") )
	#ggboxplot(test, x = "gene", y = "sces",
	#		  color = "black",  add = "jitter",width=0.25, add.params = list(color = "gene"),
	#		  xlab="Gene",ylab="SCEs/library") +
	#	stat_compare_means(comparisons = my_comparisons,label.y = c(31, 34, 37)) +
	#	stat_compare_means(label = "p.signif", method = "t.test",
	#					   ref.group = "WT")



	breakpointSummary$width=breakpointSummary$end- breakpointSummary$start

	suppressMessages(suppressWarnings(ggplot(breakpointSummary) + geom_smooth(aes(Reads_per_Mb, width,color=gene),se=F) +
									  	scale_y_log10() +
									  	theme_classic() +
									  	theme(text=element_text(size=15))+
									  	geom_hline(yintercept=10000, linetype="dashed", color = "red") +
									  	ggsave("Output/Plots/resolutionVsDepth.png")))


	sce_summary=data.frame(gene=character(),SCE=numeric(),mean_resolution=numeric(),median_resolution=numeric())
	b<- filter(breakpointSummary, gene=="BLM")
	sce_summary<- add_row(sce_summary,gene="BLM",SCE=nrow(b),mean_resolution=mean(b$width),median_resolution=median(b$width))
	r<- filter(breakpointSummary, gene=="RECQL5")
	sce_summary<- add_row(sce_summary,gene="RECQL5",SCE=nrow(r),mean_resolution=mean(r$width),median_resolution=median(r$width))
	br<- filter(breakpointSummary, gene=="BLM/RECQL5")
	sce_summary<- add_row(sce_summary,gene="BLM/RECQL5",SCE=nrow(br),mean_resolution=mean(br$width),median_resolution=median(br$width))
	w<- filter(breakpointSummary, gene=="WT")
	sce_summary<- add_row(sce_summary,gene="WT",SCE=nrow(w),mean_resolution=mean(w$width),median_resolution=median(w$width))
	write.table(sce_summary,"Output/Tables/SCE_summary.txt",quote=F,row.names = F,col.names = T,sep="\t")




	#plot here
	suppressMessages(suppressWarnings(ggplot(breakpointSummary) + stat_ecdf(aes(width,color=gene)) +
									  	scale_x_log10() +
									  	theme_classic() +
									  	ylab("SCEs Mapped (%)") +
									  	xlab("Resolution") +
									  	annotation_logticks(sides = "b") +
									  	theme(text = element_text(size=15))+
									  	geom_density(aes(width),size=1.1)+
									  	geom_vline(xintercept=median(breakpointSummary$width), linetype="dashed", color = "red") +
									  	geom_text(aes(x=5000, label=paste0("Median\n",median(width)," bp"), y=0.8))  +
									  	ggsave("Output/Plots/breakpointResolution.png")))
}
