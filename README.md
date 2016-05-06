# CBB752_Final_Project_1.3
A tool that trims reads based on quality score from a FastQ file.  Requires FastQ input file and quality score conversion table. Optionally accepts user-specified quality score trim threshold, sliding window size, minimum read length, output file type, and output file name. If no input is specified for optional parameters, the following defaults are used:   
* Threshold score: 20  
* Window size: 4  
* Minimum read length: none  
* Output file type: FastQ  
* Output file name: "output"  


## To use the python file: 
  Download "QualityTrim.py" and "qscores.txt", and optionally the example input file provided.
  
  Usage: 
  
```
python QualityTrim.py -i <input fastq file> -s <qscore file> -t <threshold cutoff score> -w <sliding window size> -m <minimum read length> -f <output file type {fastq, fasta, both}> -o <output file name>
```
  
  Example: 
  
```
python QualityTrim.py -i input.fastq -s qscores.txt -f both -o python_example_output
```

  
  Example output for the first 100 lines (25 sequences) of the example input file is included as "python_example_output.fastq" and "python_example_output.fasta"
  
  
## To use the R file: 

Download "QualityTrim.R" and "qscores.txt", and optionally the example input file provided.
  
  Usage: 
  
```Rscript QualityTrim.R -i <input fastq file> -s <qscore file> -t <threshold cutoff score> -f <output file format {fastq, fasta, both}> -m <minimum read length> -o <output file base name>```
  
  
  Example: 
  
```Rscript QualityTrim.R -i input.fastq -s qscore.txt -t 25 -f fastq -m 10 -o outputfile```
  
  An example output is shown as output_R.fasta and output_R.fastq