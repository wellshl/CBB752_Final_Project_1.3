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
# Updated.date  :  
# Updated.by    : 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Usage:      Rscript QualityTrim.R /path/to/inputfile.fastq /path/to/outputfile.fastq
# Example:    example line
# Note:       Input: fastq file, which is a text file with 1st row = header
#             
#           
#             Output: a fastq with 3' ends truncated by quality score input
#
#             example.fastq: 


# DS start of code




# [[Notes from Mtg]]
#
# Pre-paired end assembly
# 
# sliding window average for q score, ILMN
# 
# trim 3' end if average is lower
# 
# read length cutoff
# 
# default q-score but user input option
# 
# user specified output as fastq or fasta or both

# read in data file
# 
# <<< send out data file to collabs
# 
# Q scores converted to #s
# 
# # set sliding window
# 
# start sliding at 5' end
# if above threshold move windown
# else cut
# 
# score less than threshold 
# cut when ave Q below