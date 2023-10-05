import os
import subprocess

def index_ref(ref_file, output_BAM):
    # Indexes reference file using bwa index command
    cmd = ['bwa', 'index', '-p', os.path.join(output_BAM, 'ref_index'), ref_file] # Command to index ref file
    subprocess.run(cmd, check=True) # Performs indexing of ref file

def align(ref_file, fastq_preproc_file, output_BAM):
    # Performs alignment of fastq preprocessed file to ref file using bwa mem command
    sample_name = os.path.splitext(os.path.basename(fastq_preproc_file))[0] # Gets name from preprocessed FASTQ filename
    sam_file = os.path.join(output_BAM, f'{sample_name}.sam') # Defines output SAM file path
    cmd = ['bwa', 'mem', os.path.join(output_BAM, 'ref_index'), fastq_preproc_file]   # Command for alignment
    with open(sam_file, 'w') as outfile:  # Performs alignment and writes output to SAM file
        subprocess.run(cmd, stdout=outfile, check=True)


def main():
    # Specifies paths to reference and fastq files
    ref_file = 'dgorgon_reference.fa'
    fastq_dir = 'fastqs'
    output_BAM = 'BAM'
    
    if not os.path.exists('BAM'): # Creates output directory 'BAM' if doesn't exist
        os.makedirs('BAM')

    # Indexes reference file
    index_ref(ref_file, output_BAM)

    # Performs alignment for each FASTQ preprocessed file to indexed ref file
    for fastq_preproc_file in os.listdir(fastq_dir):
        if fastq_preproc_file.endswith('_trimmed.fastq'):
            fastq_preproc_file_path = os.path.join(fastq_dir, fastq_preproc_file)
            align(ref_file, fastq_preproc_file_path, output_BAM)

if __name__ == "__main__":
    main()

