# CBB752_Final_Project_1.3
A tool that trims reads based on quality score from a FastQ file.  Requires FastQ input file and quality score conversion table. Accepts user-specified quality score trim threshold, output file type, and minimum read length. If no threshold is specified, the default of 20 is used. If no output file type is specified, the default of FastQ is returned. If no minimum read length is specified, all reads are returned.


# To use the python file: 
  Download "QualityTrim.py" and "qscores.txt", and optionally the example input file provided.
  
  Usage: python QualityTrim.py -i <input fastq file> -s <qscore file> -t <threshold cutoff score> -o <output file type {fastq, fasta, both}> -m <minimum read length>
  
  To store the ouptut file: use >outputfile.txt 
  
  Example: python QualityTrim.py -i input.txt -s qscore.txt -t 25 -o fastq -m 10 > outputfile.txt
  
  Example output for the first 100 lines (25 sequences) of the example input file is shown in the folder "HLW" with and without a minimum read length of 10.
  
  
