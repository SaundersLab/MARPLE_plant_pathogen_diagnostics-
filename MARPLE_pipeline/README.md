# MARPLE plant pathogen diagnostics
Scripts for realtime surveillance of plant pathogens using the MinION sequencer (Radhakrishnan et al.,)


## Prerequisites
``` 
python 3
biopython 1.65
cufflinks 2.2.1
seqkit 0.9.1
EMBOSS 6.5.7
samtools 1.7
bwa 0.7.12
porechop 0.2.3
```

## Installation
The easiest way to set up your machine (specifications as recommended by Oxford Nanopore Technologies Ltd. for running the MinION) with the appropriate prerequisites for using the MARPLE pipeline is to use the Anaconda distribution: https://www.anaconda.com/distribution/#download-section

Once you have downloaded and installed the Anaconda distribution (instructions on the Anaconda website: https://docs.anaconda.com/anaconda/), you can install all the prerequisite software by running the following commands using the terminal:
```
conda install -c conda-forge biopython=1.65 porechop=0.2.3
conda install -c bioconda cufflinks=2.2.1 seqkit=0.9.1 emboss=6.5.7 samtools=1.7 bwa=0.7.12 
```

## Running the pipeline
With the above prerequisites installed, and this repository cloned into the folder where you want to run the analysis or by adding where the repository has been cloned into your PATH variable, run the marple_pipeline.sh script in the following format (using basecalled fastq file for a single sample as input): 

``
bash marple_pipeline.sh $fastq_file
``
