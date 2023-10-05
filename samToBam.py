import os
import subprocess

# Uses samtools view to convert SAM file to BAM file
def convert(sam_file, bam_output_dir):
    bam_file = os.path.join(bam_output_dir, os.path.basename(sam_file).replace('.sam', '.bam')) # Constructs path for output BAM file from input SAM file
    cmd = ['samtools', 'view', '-bS', sam_file, '-o', bam_file] # Command to convert SAM to BAM file
    subprocess.run(cmd, check=True) # Performs SAM to BAM file conversion

# Uses samtools sort to sort BAM file
def sort(bam_file, bam_output_dir):
    bam_filename = os.path.basename(bam_file)  # Extracts base filename from input BAM file path
    sorted_bam_filename = bam_filename.replace('.bam', '.sorted.bam')  # Generates filename for sorted BAM file
    sorted_bam_file = os.path.join(bam_output_dir, sorted_bam_filename)  # Constructs full path for sorted BAM file
    cmd = ['samtools', 'sort', '-m', '100M', '-o', sorted_bam_file, bam_file]  # Command for samtools sort 
    subprocess.run(cmd, check=True)  # Performs samtools sort command to sort BAM file

# Uses samtools index to index sorted BAM file
def index(sorted_bam_file):    
    cmd = ['samtools', 'index', sorted_bam_file]  # Command to create an index for sorted BAM file
    subprocess.run(cmd, check=True)  # Execute samtools index command to make BAM index

# Uses OS to delete files
def delete(file_list):
    for file in file_list:
        os.remove(file)

def main():
    # Specifies path to SAM files and output directory
    sam_dir = 'BAM'
    bam_output_dir = 'BAM'

    # Converts SAM files to BAM, then sorts and indexes BAM files
    for sam_file in os.listdir(sam_dir): 
        if sam_file.endswith('.sam'): # Check if the file is a SAM file
            sam_file_path = os.path.join(sam_dir, sam_file) # Constructs full path of SAM file
            convert(sam_file_path, bam_output_dir) # Converts SAM file to BAM format
            sort(os.path.join(bam_output_dir, os.path.basename(sam_file).replace('.sam', '.bam')), bam_output_dir) # Sort BAM files 
            index(os.path.join(bam_output_dir, os.path.basename(sam_file).replace('.sam', '.sorted.bam')))  # Create index for sorted BAM files

    # Get list of all .sam files and unsorted .bam files in BAM directory
    sam_to_delete = [os.path.join(sam_dir, file) for file in os.listdir(sam_dir) if file.endswith('.sam')]
    bam_to_delete = [os.path.join(bam_output_dir, file) for file in os.listdir(bam_output_dir) if file.endswith('.bam') and not file.endswith('.sorted.bam')]

    # Deletes .sam and unsorted .bam files
    delete(sam_to_delete)
    delete(bam_to_delete)

if __name__ == "__main__":
    main()

