### Title: QualityTrim.py
### Version: 1.0
### Purpose: A tool that trims the 3' end of reads using sliding window average Q-score threshold
### Created Date: 01 May 2016
### Created By: Heather Wells

### Usage: python QualityTrim.py -i <input fastq file> -s <qscore file> -t <threshold cutoff score> -m <minimum read length> -f <output file type {fastq, fasta, both}> -o <output file name>
### Example: python QualityTrim.py -i input.txt -s qscore.txt -t 25 -o fastq -m 10 > outputfile.txt

### Notes: Requires input file in fastq format and score file in "|" separated format.
###        If no threshold cutoff score is specified the default of 40 is used.
###        If no output file type is specified the default fastq is used.
###        If no minimum read length is specified, no minimum is used.
###        To save output files, use >outputfile.txt 

#import arguments
import argparse
parser=argparse.ArgumentParser(description='QualityTrim')
parser.add_argument('-i', '--input', help='input file name', required=True)
parser.add_argument('-s', '--score', help='score file name', required=True)
parser.add_argument('-t', '--threshold', help='threshold cutoff score', default=20, type=int)
parser.add_argument('-w', '--window', help='window size', default=4, type=int)
parser.add_argument('-m', '--minimum', help='minimum read length', default=0, type=int)
parser.add_argument('-f', '--outputtype', help='output file type, must be one of the following: {fastq, fasta, both}', default='fastq')
parser.add_argument('-o', '--outputname', help='output file name', default="output")
args=parser.parse_args()

import re
import numpy


def Trim(inputfile,scorefile,threshold,outputtype):

    #open and save inputs
    with open(args.input) as input:
        input=input.read().splitlines()[0:100]
    with open(args.score) as scores:
        scores=numpy.loadtxt(scores,delimiter="|",dtype=str,comments=None)

    #initialize empty output lists
    fastq=[]
    fasta=[]

    #for loop with step size of 4 to analyze one sequence chunk at a time
    for i in range(0,len(input),4):
        
        #the ascii qscore is the 4th (i+3) of each chunk of 4 lines
        q_ascii_i=list(input[i+3])

        #convert ascii code into numeric scores
        #for each character, find the match in the scores array and add to q_score_i list
        #then, for each q_score_i list, average the scores 4 at a time (sliding window)
        #if the average is less than the specified or default (20) threshold, break the loop and store the cut value as the window midpoint
        #if the average is greater than the default for the entire loop, set cut to length of q_score_i (keeps entire sequence)
        #the sequence is the 2nd (i+1) of each chunk of 4 lines, only keep sequence up to the cut point
        q_score_i=[]
        window=args.window
        for each in q_ascii_i:
            q_score_i.extend(map(int,scores[numpy.where(scores==each)[0],1]))
        for j in range(len(q_score_i)-window):
            avg=sum(q_score_i[j:j+window])/window
            if avg<args.threshold:
                cut=j+window/2
                break
            else: cut=len(q_score_i)
        seq_i=input[i+1][0:cut]
        
        #extend fasta and fastq file types with correct format only if length minimum is met (defaults with no min)
        if len(seq_i)>args.minimum:
            fastq.extend(input[i]+"\n"+seq_i+"\n+\n"+input[i+3]+"\n")
            fasta.extend(">"+input[i]+"\n"+seq_i+"\n")

    fastq="".join(fastq)
    fasta="".join(fasta)

    #save user-specified or default (fastq) output
    if args.outputtype=='fastq':
        fastq_file=open(args.outputname+".fastq","w")
        fastq_file.write(fastq)
        fastq_file.close()
    if args.outputtype=='fasta':
        fasta_file=open(args.outputname+".fasta","w")
        fasta_file.write(fasta)
        fasta_file.close()
    if args.outputtype=='both':
        fastq_file=open(args.outputname+".fastq","w")
        fastq_file.write(fastq)
        fastq_file.close()
        fasta_file=open(args.outputname+".fasta","w")
        fasta_file.write(fasta)
        fasta_file.close()

Trim(args.input,args.score,args.threshold,args.outputtype)
