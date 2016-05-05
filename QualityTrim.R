#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title :  QualityTrim.R 
# Version : 1.0
#
# Purpose : A tool that trims the 3' end of reads using sliding window average Q-score threshold
#  
# Version Notes : 
#
# Created.date  : 27 Apr 2016
# Created.by    : Dan Spakowicz
# Updated.date  : 05 May 2016 
# Updated.by    : DS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Usage: Rscript QualityTrim.py -i <input fastq file> -s <qscore file> -t <threshold cutoff score> -o <output file type {fastq, fasta, both}> -m <minimum read length>
### Example: Rscript QualityTrim.py -i input.txt -s qscore.txt -t 25 -o fastq -m 10 > outputfile.txt

### Notes: Requires input file in fastq format and score file in "|" separated format.
###        If no threshold cutoff score is specified the default of 40 is used.
###        If no output file type is specified the default fastq is used.
###        If no minimum read length is specified, no minimum is used.
###        To save output files, use > outputfile.txt 


rm(list=ls())
oldw <- getOption("warn")
options(warn = -1)

# Load the required packages
list.of.packages <- c("optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(optparse)

# set arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input fastq file name", metavar="character"),
  make_option(c("-s", "--score"), type="character", default="./qscores.txt", 
              help="score file name", metavar="character"),
  make_option(c("-t", "--threshold"), type="character", default=20, 
              help="threshold cutoff score", metavar="character"),
  make_option(c("-m", "--minimum"), type="character", default=0, 
              help="minimum read length", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              default="fastq", 
              help="output file type, must be one of the following: {fastq, fasta, both}, [default= %default]", 
              metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# infile <- "example_100.fastq"
# scorefile <- "qscores.txt"
# threshold <- 20
# minimum <- 20
# outputtype <- "fastq"

# Trim function
Trim <- function(infile, scorefile, threshold, minimum, outputtype) {
  
  # read in and format fastq file
  input <- readLines(infile)
  fastq.split <- strsplit(as.character(input), split = "")
  
  # read in the quality scores file
  scorefile <- read.table(scorefile, sep = "|", quote = "", comment.char = "")
  
  # list of just the qual score lines
  qualLines <- seq(4, length(input), by = 4)
  qualscores <- fastq.split[qualLines]
  
  # convert ascii qual scores to phred scores
  phred.scores <- vector("list", length(qualscores))
  for ( i in 1:length(qualscores)) {
    for (j in 1:length(qualscores[[i]])){
      phred.scores[[i]][j] <- scorefile$V2[which(qualscores[[i]][j]==scorefile$V1)]
    }
  }
  
  # sliding window average over the phred scores
  cut.logical <- vector("list", length(phred.scores))
  for ( i in 1:length(phred.scores)) {
    for (j in 1:length(phred.scores[[i]])){
      cut.logical[[i]][j] <- sum(phred.scores[[i]][j:(j+3)])/4 < threshold
        }
    }
  
  # find the location of the first average less than the threshold (TRUE)
  cut.loc <- vector("list", length(cut.logical))
  for ( i in 1:length(cut.logical)) {
    cut.loc[[i]] <- min(which(cut.logical[[i]] == TRUE))
    if(cut.loc[[i]] == Inf){cut.loc[[i]] = length(cut.logical[[i]])}
  }
 
  # trim sequence and qual scores
  cut.loc.vec <- unlist(cut.loc)
  for (i in 1:length(cut.loc.vec)){
    fastq.split[[(i*4-2)]] <- fastq.split[[(i*4-2)]][1:(cut.loc.vec[i]+1)]
    fastq.split[[(i*4)]] <- fastq.split[[(i*4)]][1:(cut.loc.vec[i]+1)]
  }
  
  # remove sequences shorter than the minimum length
  minlength <- vector("list", length(qualLines))
  for (i in qualLines){
    minlength[[i]] <- length(fastq.split[[i]]) > minimum
  }
  for (i in qualLines){
    minlength[[(i-3)]] <- minlength[[i]]
    minlength[[(i-2)]] <- minlength[[i]]
    minlength[[(i-1)]] <- minlength[[i]]
  }
  minlength.vec <- unlist(minlength)
  fastq.split.min <- fastq.split[minlength.vec]
  
  # output as fastq, fasta or both
  if ( outputtype == "fasta" ){
    
    a <- seq(from = 1, to = length(fastq.split.min), by = 4)
    b <- seq(from = 2, to = length(fastq.split.min), by = 4)               
    fastaLines <- c( matrix(c(a,b), nrow=2, byrow=TRUE) )
    fasta <- fastq.split.min[fastaLines]
    for ( i in 1: length(fasta)){
      if(fasta[[i]][1]=="@"){
        fasta[[i]][1] <- ">"
      }
    }
    fileFasta<-file("output.fasta")
    writeLines(unlist(lapply(fasta, paste, collapse="")), fileFasta)
    close(fileFasta)
  }
  
  if ( outputtype == "fastq" ){
    fileFastq<-file("output_R.fastq")
    writeLines(unlist(lapply(fastq.split.min, paste, collapse="")), fileFastq)
    close(fileFastq)
  }
  if ( outputtype == "both" ){
    fileFastq<-file("output_R.fastq")
    writeLines(unlist(lapply(fastq.split.min, paste, collapse="")), fileFastq)
    close(fileFastq)
    a <- seq(from = 1, to = length(fastq.split.min), by = 4)
    b <- seq(from = 2, to = length(fastq.split.min), by = 4)               
    fastaLines <- c( matrix(c(a,b), nrow=2, byrow=TRUE) )
    fasta <- fastq.split.min[fastaLines]
    for ( i in 1: length(fasta)){
      if(fasta[[i]][1]=="@"){
        fasta[[i]][1] <- ">"
      }
    }
    fileFasta<-file("output_R.fasta")
    writeLines(unlist(lapply(fasta, paste, collapse="")), fileFasta)
    close(fileFasta)
  }
}

Trim(opt$input, opt$score, opt$threshold, opt$minimum, opt$out)

options(warn = oldw)  
