Title: Bioinformatics pipeline for variant analysis 
===================================================
This repository comprises a collection of Python programs designed to automate a bioinformatics pipeline for variant analysis. These scripts have been validated to run on Linux systems to process DNA sequence data, align it to a reference genome, identify genetic variants, and produce a summary report.
Purpose  
-------
preprocessFastq.py:
The preprocessFastq.py script preprocesses pooled raw FASTQ files containing DNA sequences. It performs tasks such as demultiplexing pooled samples, trimming adapter barcodes, and filtering sequences based on quality. The script uses ParseFastQ  to read and manipulate FASTQ sequencing data. It identifies sample-specific barcodes to separate sequences into individual FASTQ files. It uses filtering thresholds to retain only high-quality, non-degraded sequences.

alignFastq.py:
The alignFastq.py script aligns the preprocessed DNA sequences from the FASTQ files to the reference genome. It uses BWA (Burrows-Wheeler Aligner) to align sequences to the reference genome. The script generates SAM (Sequence Alignment/Map) files that contain the data for the alignment of each DNA sequence to the reference.

samToBam.py:
The samToBam.py script converts the SAM files produced by alignFastq.py into the BAM (Binary Alignment/Map) file. BAM files are more efficient than SAM files because they reduce storage requirements and allow for quicker data access. The script uses the pysam library to manipulate SAM and BAM files, for file conversion, sorting, and indexing.

findVariant.py:
The findVariant.py script analyzes aligned BAM files to identify genetic variants. It scans through the aligned reads, compares them to the reference genome, and detects base positions with variants. It calculates the frequency of each variant.

pipeline.py:
The pipeline.py script serves as the coordinator for the entire analysis pipeline. It executes the preceding scripts in a sequential manner, passing the required inputs and outputs. The script automates the entire variant analysis pipeline. It generates a report, report.txt, summarizing the mutations and frequencies in each sample.

Execute the scripts using Terminal in Ubuntu/Linux:  
---------------------------------------------------
1. Check that scripts are in the correct directory.

The preprocessFastq.py, alignFastq.py, findVariant.py, samToBam.py, and pipeline.py scripts should be located in the working directory.

The dgorgon_reference.fa, hawkins_pooled_sequences.fastq, and harrington_clinical_data.txt files should be located in the working directory.

2. Install dependencies

pip install pysam
sudo apt install bwa
sudo apt install samtools

3. Using Terminal, navigate to the working directory.  
 
4. Give the scripts executable permissions using the following commands:  

chmod +x preprocessFastq.py
chmod +x alignFastq.py 
chmod +x samToBam.py
chmod +x findVariant.py
chmod +x pipeline.py
 
5. Execute the script with the following command:  
 
python3 pipeline.py -f hawkins_pooled_sequences.fastq 
 
6. The following outputs are generated:  
 
A folder named `fastqs` containing demultiplexed fastqs for each sample. 

A folder named `BAMS` containing sorted.bam and .bai files for each sample. 

A text file named report.txt containing sample names, mold colors, mutations, mutation locations, and mutation frequencies.

 
