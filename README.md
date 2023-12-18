
# RUSH - regulation uncovered by steric hindrance

Data analysis on gRNA sequencing results mainly uses two python scripts, parse.py and rush.py.

parse.py processes gRNA amplicon sequencing library to obtain individual gRNA counts for each screen, which includes the following steps:

    1. concatenate and unzip fastq.gz files (.fastq.gz -> .fastq)
    2. trim Illumina sequences (.fastq -> trimmed.fastq)
    3. remove primer adapters (trimmed.fastq -> ready.fastq)
    4. count perfect gRNAs (ready.fastq -> library_count.csv)

rush.py handles gRNA count csv files from all the screens to generate results in summary csv files.


Relevant file descriptions are listed below:

    count_spacers_py3.py - module that finds gRNA sequences (parse.py)
    TruSeq3-SE.fa - fasta file needed for trimmomatic that trims adapter sequences (parse.py)
    library_sequences.csv - csv file containing gRNA library sequences (parse.py)

    utils.py - module containing functions for calculations (rush.py)
    oligo_library_sequences.xlsx - excel file containing information about the gRNA library (rush.py)

------------------------------------------------------------------------------------------------------
Several packages need to be installed before running the scripts - 

    conda install -c bioconda trimmomatic
    
    # cupadapt 3.4 with Python 3.8.13
    conda install -c bioconda cutadapt
    conda install -c conda-forge biopython


The scripts requires specific data organization. For parse.py, relevant data files (TruSeq3-SE.fa, library_sequences.csv, fastq.gz files) and scripts (parse.py, count_spacers_py3.py) are stored under the same directory (homeDir).

    homeDir
        -count_spacers_py3.py
            > unsorted_library_count.csv
            > top5_library_count.csv
            > bot5_library_count.csv
        -TruSeq3-SE.fa
        -library_sequences.csv
        -unsorted101122_S1_L001_R1_001.fastq.gz
        -unsorted101122_S1_L002_R1_001.fastq.gz
        -unsorted101122_S1_L003_R1_001.fastq.gz
        -unsorted101122_S1_L004_R1_001.fastq.gz
        -top5...
        -bot5...
        
For rush.py, each screen folder contains three count csv files (unsorted, top5, bot5). In our case, four different spliceRUSH screens correspond to four different folders and twelve csv files. The folder name should contain "dualEx" or "dualIn" (not case-sensitive) for the script to handle the files properly.

    inputDir
        -DualEX-SMN2-Rep1-101122
            > unsorted_library_count.csv
            > top5_library_count.csv
            > bot5_library_count.csv
        -DualEX-SMN2-Rep2-102122
        -DualIN-SMN2-Rep1-091721
        -DualIN-SMN2-Rep2-020722
        -oligo_library_sequences.xlsx
      
------------------------------------------------------------------------------------------------------
Code example:

    python parse.py -d homeDir --libcsv library_sequences.csv -u unsorted -t top5 -b bot5
    python rush.py -i inputDir -l oligo_library_sequences.xlsx

