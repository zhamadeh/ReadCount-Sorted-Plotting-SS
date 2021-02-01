

summarizeBreaks <- function(breakpoints) {
    breaks <- breakpoints$breaks
    confint <- breakpoints$confint
    if (length(breaks) > 0 & length(confint) == length(breaks)) {
        breaks.df <- as(breaks, 'data.frame')
        confint.df <- as(confint, 'data.frame')
        breaks.df <- breaks.df[,c('seqnames','start','end')]
        confint.df <- confint.df[,c('start','end','genoT',"deltaW")]
        breaksSummary <- cbind(breaks.df, confint.df)
        names(breaksSummary) <- c('seqnames','start','end','CI.start','CI.end','genoT',"deltaW")
        return(breaksSummary)
    } else if (length(breaks) > 0 & length(confint) == 0) {
        breaks.df <- as(breaks, 'data.frame')
        breaksSummary <- breaks.df[,c('seqnames','start','end','genoT')]
        return(breaksSummary)
    } else {
        return(NULL)
    }
}
