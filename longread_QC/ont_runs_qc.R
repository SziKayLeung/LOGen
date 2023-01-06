# Szi Kay Leung
# Functions adapted from https://github.com/nanoporetech/ont_tutorial_basicqc.git Rmarkdown

suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(extrafont))


# prepare_summary_file <sequencing_summary.txt>
# Read in txt file and extract the number of passed and failed sequences
prepare_summary_file <- function(inputFile){
  # Using fread for fast and friendly import of sequence_summary file
  # no definition of column types to speed import and allow for merged files
  # could be worthwhile to fread(select=) to select only a subset of columns - this could
  # preclude e.g. barcode data or different versions?
  sequencedata <- data.table::fread(inputFile, stringsAsFactors=FALSE)
  # remove the redundant headers from merged files
  if (length(which(sequencedata[,1]=="filename")) > 0) {
    sequencedata <- sequencedata[-which(sequencedata[,1]=="filename"),]
  }
  # coerce the columns used in analytics into more appropriate data-types 
  sequencedata$channel<-as.numeric(sequencedata$channel)
  sequencedata$start_time<-as.numeric(sequencedata$start_time)
  sequencedata$duration<-as.numeric(sequencedata$duration)
  sequencedata$num_events<-as.numeric(sequencedata$num_events)
  sequencedata$sequence_length_template<-as.numeric(sequencedata$sequence_length_template)
  sequencedata$mean_qscore_template<-as.numeric(sequencedata$mean_qscore_template)
  # passes_filtering is a useful flag; but there are examples of sequencing_summary.txt where this 
  # is not present - https://github.com/a-slide/pycoQC/blob/master/pycoQC/data/sequencing_summary_1D_DNA_Albacore_1.2.1.txt
  if (! "passes_filtering" %in% colnames(sequencedata)) {
    # set all of the reads to pass? apply a cutoff?
    sequencedata$passes_filtering <- TRUE
  } else {
    sequencedata$passes_filtering <- as.logical(sequencedata$passes_filtering)
  }
  # create a convenient separation of pass and fail ...
  passedSeqs <- sequencedata[which(sequencedata$passes_filtering), ]
  failedSeqs <- sequencedata[which(!sequencedata$passes_filtering), ]
  
  output <- list(sequencedata, passedSeqs, failedSeqs)
  names(output) <- c("sequencedata","passedSeqs","failedSeqs")
  return(output)
}

key_metrics <- function(df, passedSeqs){
  readCount <- formatC(nrow(df$sequencedata), big.mark=",")
  totalBases = sum(df$sequencedata$sequence_length_template,na.rm=T)/10^9
  passedbasecalled <- formatC(nrow(passedSeqs), big.mark=",")
  perc_passedbasecalled <- round(nrow(passedSeqs) / nrow(df$sequencedata) * 100, 1)
  passedBases <- round(sum(df$passedSeqs$sequence_length_template,na.rm=T)/10^9,2)
  perc_passedBases <- round(passedBases / totalBases * 100, 1)
  gigabases <- round(totalBases,2)
  
  metric.data <- data.frame(readCount,totalBases, gigabases,passedbasecalled, perc_passedbasecalled, passedBases, perc_passedBases)
  metric.data <- gather(metric.data, metric,value, readCount:perc_passedBases, factor_key=TRUE)
  
  metric.data$metric <- c("Total Number of Reads Basecalled","Total number of Bases","GB","Number of Reads that passed Basecalling","Percentage of Number of Reads that passed Basecalling","Total Number of Gb of DNA sequence in passed reads","Percentage of Gb passed") 
  
  return(metric.data)
}

channel_create <- function(sequencedata){
  # create an empty read count container ... MinION or PromethION??
  # https://gist.github.com/roblanf/df47b9748c3aae00809cc675aca79989
  # build the map for R9.5 flowcell, as a long-form dataframe that translates
  # channels into rows and columns on the flowcell. Good for plotting in R.
  p1 = data.frame(channel=33:64, row=rep(1:4, each=8), col=rep(1:8, 4))
  p2 = data.frame(channel=481:512, row=rep(5:8, each=8), col=rep(1:8, 4))
  p3 = data.frame(channel=417:448, row=rep(9:12, each=8), col=rep(1:8, 4))
  p4 = data.frame(channel=353:384, row=rep(13:16, each=8), col=rep(1:8, 4))
  p5 = data.frame(channel=289:320, row=rep(17:20, each=8), col=rep(1:8, 4))
  p6 = data.frame(channel=225:256, row=rep(21:24, each=8), col=rep(1:8, 4))
  p7 = data.frame(channel=161:192, row=rep(25:28, each=8), col=rep(1:8, 4))
  p8 = data.frame(channel=97:128, row=rep(29:32, each=8), col=rep(1:8, 4))
  q1 = data.frame(channel=1:32, row=rep(1:4, each=8), col=rep(16:9, 4))
  q2 = data.frame(channel=449:480, row=rep(5:8, each=8), col=rep(16:9, 4))
  q3 = data.frame(channel=385:416, row=rep(9:12, each=8), col=rep(16:9, 4))
  q4 = data.frame(channel=321:352, row=rep(13:16, each=8), col=rep(16:9, 4))
  q5 = data.frame(channel=257:288, row=rep(17:20, each=8), col=rep(16:9, 4))
  q6 = data.frame(channel=193:224, row=rep(21:24, each=8), col=rep(16:9, 4))
  q7 = data.frame(channel=129:160, row=rep(25:28, each=8), col=rep(16:9, 4))
  q8 = data.frame(channel=65:96, row=rep(29:32, each=8), col=rep(16:9, 4))
  # long form as a data frame, i.e. map$channel[[1]] returns 33
  channelMap = rbind(p1, p2, p3, p4, p5, p6, p7, p8, q1, q2, q3, q4, q5, q6, q7, q8)
  hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab') #RdPu, Oranges, Greens, YlOrRd, Purples
  channelCounts <- as.data.frame(matrix(rep(0, 512), ncol=1))
  channelCountRaw <- as.data.frame(table(unlist(sequencedata[, "channel"])), row.names=1)
  channelCounts[row.names(channelCountRaw),] <- channelCountRaw[,1]
  #channelMap <- cbind(channelMap[channelMap$channel,], frequency=channelCounts[channelMap$channel,])
  channelMap <- merge(channelMap, channelCounts, by.x="channel", by.y=0)
  colnames(channelMap)[4]<-"count"
  channelMapMatrix <- reshape2::acast(channelMap, col ~ row, value.var = "count")
  
  p <- ggplot(channelMap, aes(x = row, y = col, fill = count)) +
    geom_tile() +
    geom_text(data=channelMap,aes(x=row, y=col,label=count,color=count),show.legend = F, size=2.5) +
    scale_x_discrete(breaks=NULL) +
    scale_y_discrete(breaks=NULL) +
    coord_equal() +
    scale_fill_gradientn(colours = hm.palette(100), name = "Number of Sequences") +
    scale_color_gradient2(low = hm.palette(100), high = hm.palette(1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="bottom",
          legend.key.width=unit(5.6,"cm"))
  
  return(p)
  
}

quality_length <- function(df){
  lenSorted <- rev(sort(df$passedSeqs$sequence_length_template))
  passedMeanLength = round(mean(lenSorted), digits = 0)
  N50 <- lenSorted[cumsum(lenSorted) >= sum(lenSorted)*0.5][1]
  passedMeanQ = round(mean(df$passedSeqs$mean_qscore_template), digits = 1)
  failedMeanQ = round(mean(df$failedSeqs$mean_qscore_template), digits = 1)
  
  metric <-  c("Mean Read Length (nt)","N50","Mean Read Quality (QV)","Mean Failed QV","Longest Read")
  value <- c(passedMeanLength, N50, passedMeanQ, failedMeanQ, prettyNum(max(df$passedSeqs$sequence_length_template), big.mark=","))
  output <- data.frame(metric,value)
  
  return(output)
}

read_length <- function(df){
  # https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
  roundUpNice <- function(x, nice=seq(from=1, to=10, by=0.25)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  # pick a friendly upper limit to render sequence lengths into a histogram
  # here we're aiming for a robustly rounded up 97.5 quantile of the data (skip a few outliers ...)
  # an ideal histogram will have 40 or so bins
  histogramBinCount <- 40
  upperLimit <- roundUpNice(as.numeric(quantile(x=df$sequencedata$sequence_length_template, probs=c(0.975))))
  breakVal = roundUpNice(upperLimit / histogramBinCount)
  breaks <- seq(0, to=upperLimit, by=breakVal)
  binAssignments <- cut(df$sequencedata$sequence_length_template, breaks, include.lowest=TRUE, right=FALSE)
  scrapeBinnedReads <- function(level, qcpass) {
    length(subset(df$sequencedata[which(binAssignments == level), ], passes_filtering==qcpass)$sequence_length_template)
  }
  
  lenSorted <- rev(sort(df$passedSeqs$sequence_length_template))
  passedMeanLength = round(mean(lenSorted), digits = 0)
  N50 <- lenSorted[cumsum(lenSorted) >= sum(lenSorted)*0.5][1]
  
  
  passedBinnedReads <- unlist(lapply(levels(binAssignments), scrapeBinnedReads, qcpass=TRUE))
  failedBinnedReads <- unlist(lapply(levels(binAssignments), scrapeBinnedReads, qcpass=FALSE))
  binnedReadDist <- data.frame(length=head(breaks, -1), pass=passedBinnedReads, fail=failedBinnedReads)
  binnedReadMelt <- reshape2::melt(binnedReadDist, id.vars=c("length"))
  
  p <- ggplot(binnedReadMelt, aes(x=length, fill=variable, y=value)) +
    geom_bar(stat="identity") +
    xlab("Read length (kb)") + ylab("Number of Reads (Thousand)") +
    scale_fill_manual("QC", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    scale_x_continuous(labels = ks , limits=c(-breakVal,upperLimit), breaks=pretty(df$passedSeqs$sequence_length_template,n=40)) +
    scale_y_continuous(labels = ks) +
    #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #labs(title="Histogram showing distribution of read lengths across quality passing sequences", fill="QV filter")+
    geom_vline(xintercept = N50, size = 1) +
    annotate("text", x=N50, y=max(passedBinnedReads + failedBinnedReads), label = " N50", hjust=0, colour="SteelBlue") +
    geom_vline(xintercept = passedMeanLength, size = 1) +
    annotate("text", x=passedMeanLength, y=max(passedBinnedReads + failedBinnedReads), label = " Mean", hjust=0, colour="SteelBlue") + mytheme
  
  return(p)
}

QC_score <- function(df){
  p <- ggplot(df$sequencedata, aes(x=mean_qscore_template, fill=passes_filtering)) + 
    geom_histogram(breaks=seq(from=0, to=15, by=0.1)) +
    scale_fill_manual(name="QC", values=c("TRUE"=brewer.pal(6, "Paired")[2], 
                                          "FALSE"=brewer.pal(6, "Paired")[1]), labels=c( "pass", "fail"), breaks=c("TRUE", "FALSE")) +
    labs(x = "Read mean Q-score", y = "Number of Reads (Thousand)") + mytheme + scale_y_continuous(labels = ks)
  
  return(p)
}

QC_contour_plot <- function(df){
  qcThreshold <- 7
  binFilter <- 5
  # prepare the density plot, but do not render
  lq_dens <- ggplot(df$sequencedata, aes(log10(sequence_length_template), mean_qscore_template)) + geom_bin2d(bins=100)
  # extract the density map from the plot
  lq_dens_counts <- ggplot_build(lq_dens)$data[[1]]
  if (binFilter > 0) {
    # remove the bins from the density map that do not contain sequence count above threshold 
    lq_dens_counts <- lq_dens_counts[-which(lq_dens_counts$count <= binFilter),]
  }
  # directly plot this modified density map (stat=="identity")
  p <- ggplot(lq_dens_counts) + 
    geom_bin2d(aes(x,y,fill=count), stat="identity") +
    scale_fill_distiller(palette="Blues", trans="reverse") + 
    geom_hline(yintercept = qcThreshold, size = 1) + 
    scale_x_continuous(breaks = c(1,2,3,4,5), label = c("10", "100", "1000", "10,000", "100,000")) +
    #annotation_logticks(base = 10, sides = "b", scaled = TRUE) 
    labs(y = "Read Mean Q-score", x = "Read Length (bp)", fill = "Number of Reads (Thousand)") + mytheme +
    theme(legend.position = c(0.15,0.9))
  #labs(title="Contour Plot showing distribution of quality scores against log10 read lengths (all reads)")
  
  
  return(p)
}

temporal <- function(df){
  scaling <- 1
  df$sequencedata$start_time <- df$sequencedata$start_time - min(df$sequencedata$start_time)
  df$sequencedata$start_time <- df$sequencedata$start_time / scaling
  
  # assuming a 72 hour run, 5 minute intervals
  sampleHours = 48
  sampleIntervalMinutes = 60
  breaks = seq(0, sampleHours*60*60, by=60*sampleIntervalMinutes)
  binass <- findInterval(df$sequencedata$start_time, breaks)
  mergeItPerHour <- function(interval, binnedAssignments, filter) {
    totalbases = 0
    if (length(which(binnedAssignments==interval))>0) {
      subset <- df$sequencedata[which(binnedAssignments==interval), ]
      if (length(which(subset$passes_filtering == filter)) > 0) {
        totalbases = sum(subset[which(subset$passes_filtering == filter), "sequence_length_template"])
      }
    }
    # need to scale what is being returned - totalbases value is total bases within an interval (sampleIntervalMinutes)
    return(totalbases / 1e9 / sampleIntervalMinutes * 60)
  }
  binnedTemporalDataPerHour <- data.frame(
    cbind(
      time=breaks,
      pass=unlist(lapply(seq(breaks), mergeItPerHour, binnedAssignments=binass,filter=TRUE)),
      fail=unlist(lapply(seq(breaks), mergeItPerHour, binnedAssignments=binass, filter=FALSE))
    )
  )
  binnedTemporalDataPerHour$time <- binnedTemporalDataPerHour$time / 60 / 60
  
  return(binnedTemporalDataPerHour)
}

temporal_plot <- function(binnedTemporalDataPerHour){
  p <- ggplot(binnedTemporalDataPerHour, aes(time)) +
    geom_line(aes(y = fail, colour = "fail"), size=1) + 
    geom_line(aes(y = pass, colour = "pass"), size=1) +
    scale_color_manual(name="QV", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    xlab("Time (hours)") + 
    ylab("Gigabases sequenced per hour") + mytheme
  #labs(title="Plot showing sequence throughput against time")
  
  return(p)
  
}
temporal_cum <- function(binnedTemporalDataPerHour,df){
  # assuming a 48 hour run, 5 minute intervals
  sampleHours = 48
  sampleIntervalMinutes = 60
  # binnedTemporalDataPerHour is scaled to Gbp per hour - rescale to raw for cumulative plotting
  binnedTemporalDataPerHour$pass <- binnedTemporalDataPerHour$pass / 60 * sampleIntervalMinutes
  binnedTemporalDataPerHour$fail <- binnedTemporalDataPerHour$fail / 60 * sampleIntervalMinutes
  # https://stackoverflow.com/questions/31404679/can-ggplot2-find-the-intersections-or-is-there-any-other-neat-way
  acquireTimePoints <- which(binnedTemporalDataPerHour$pass > 0)
  targetInterpolate <- approxfun(x=binnedTemporalDataPerHour[acquireTimePoints, "time"], y=cumsum(binnedTemporalDataPerHour[acquireTimePoints, "pass"]))
  base50 <- sum(df$passedSeqs$sequence_length_template)/1e9*0.5
  base90 <- sum(df$passedSeqs$sequence_length_template)/1e9*0.9
  T50 <- optimize(function(t0) abs(targetInterpolate(t0) - base50), 
                  interval = range(binnedTemporalDataPerHour[acquireTimePoints, "time"]))
  T90 <- optimize(function(t0) abs(targetInterpolate(t0) - base90), 
                  interval = range(binnedTemporalDataPerHour[acquireTimePoints, "time"]))
  p <- ggplot(binnedTemporalDataPerHour, aes(time)) +
    geom_line(aes(y = cumsum(fail), colour = "fail"), size=1) + 
    geom_line(aes(y = cumsum(pass), colour = "pass"), size=1) +
    scale_color_manual(name="QV", values=c("fail"=brewer.pal(6, "Paired")[1], "pass"=brewer.pal(6, "Paired")[2])) +
    geom_segment(x=T50$minimum, y=0, xend=T50$minimum, yend=base50, colour="darkgray", size=1) +
    geom_segment(x=0, y=base50, xend=T50$minimum, yend=base50, colour="darkgray", size=1) +
    annotate("text", x=T50$minimum, y=base50, label=" T50", vjust=1, hjust=0, colour="SteelBlue") +
    geom_segment(x=T90$minimum, y=0, xend=T90$minimum, yend=base90, colour="darkgray", size=1) +
    geom_segment(x=0, y=base90, xend=T90$minimum, yend=base90, colour="darkgray", size=1) +
    annotate("text", x=T90$minimum, y=base90, label=" T90", vjust=1, hjust=0, colour="SteelBlue") +
    xlab("Time (hours)") + 
    ylab("Number of bases sequenced (Gb)") + mytheme
  #labs(title="Plot showing cumulative bases sequenced against time")
  
  return(p)
}

speedtime <- function(df){
  # https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
  roundUpNice <- function(x, nice=seq(from=1, to=10, by=0.25)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  # pick a friendly upper limit to render sequence lengths into a histogram
  # here we're aiming for a robustly rounded up 97.5 quantile of the data (skip a few outliers ...)
  upperLimit <- roundUpNice(as.numeric(quantile(x=df$sequencedata$sequence_length_template, probs=c(0.975))))
  # an ideal histogram will have 40 or so bins
  scaling <- 1
  histogramBinCount <- 40
  breakVal = roundUpNice(upperLimit / histogramBinCount)
  breaks <- seq(0, to=upperLimit, by=breakVal)
  
  binass <- findInterval(df$sequencedata$start_time, breaks)
  speedTime <- data.frame(segment=binass, rate=df$sequencedata$sequence_length_template / (df$sequencedata$duration/scaling))
  p <- ggplot(speedTime, aes(x=segment, y=rate, group=segment)) + geom_boxplot(fill="steelblue", outlier.shape=NA) +scale_x_continuous(name="Time (hours)") + ylab("Sequencing rate (bases per second)") + mytheme
  
  return(p)
}

AllQCPlots <- function(dat){
  
  #dat <- prepare_summary_file(sequencing_summary_input)
  time <- temporal(dat)
  
  pchannel <- channel_create(dat$sequencedata) + 
    theme(legend.position = "right",legend.key.width = unit(1, 'cm'),legend.key.height = unit(2, 'cm'))
  plength <- read_length(dat) 
  pscore <- QC_score(dat) + mytheme + theme(legend.position = "none") 
  pcontour <- QC_contour_plot(dat) + mytheme 
  ptemp <- temporal_plot(time) + mytheme
  ptemp_cum <- temporal_cum(time,dat) + mytheme 
  pspeedtime <- speedtime(dat)
  
  output <- list(pchannel,plength,pscore,pcontour,ptemp,ptemp_cum,pspeedtime)
  names(output) <- c("pchannel","plength","pscore","pcontour","ptemp","ptemp_cum","pspeedtime")
  
  return(output)
}