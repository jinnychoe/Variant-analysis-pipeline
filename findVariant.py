import pysam
import os

# Calls pileup function for BAM file to calculate nucleotide frequencies
def pileup(bam_file):
    samfile = pysam.AlignmentFile(bam_file, "rb") # Open BAM file for reading
    
    ntdict = {} # Initializes dictionary to store mutation data for each position
    total_reads = 0 # Initializes variable to count total number of reads
    
    # Loops through each pileup column (each position along reference seq) in BAM file
    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.pos  # Gets position along reference
        
        ntdict[pos] = {} # Initializes nested dictionary for this position to store base frequencies
        
        # Loops through each read in pileup column
        for pileupread in pileupcolumn.pileups:
            total_reads += 1  # Increments total read count
            
            if not pileupread.is_del and not pileupread.is_refskip: # Checks if current read is not a deletion or reference skip
                base = pileupread.alignment.query_sequence[pileupread.query_position] # Gets base at query position of current read 
                # Checks if base is already in dictionary for this position
                if base in ntdict[pos]:
                    ntdict[pos][base] += 1  # Increments count for this base
                else:
                    ntdict[pos][base] = 1  # Initializes count for this base to 1

    samfile.close() # Closes BAM file
        
    return ntdict, total_reads # Returns dictionary containing mutation data and total read count

    
# Calculates mutation frequencies in BAM folder
def calc_mut_freq(ntdict, total_reads, ref_sequence):
    mut_freq = {}  # Creates a dictionary to store mutation frequencies for each position

    # Loops through each position and its corresponding base frequencies
    for pos, bases in ntdict.items():
        total_coverage = sum(bases.values())  # Calculates total coverage at the position

        # Gets reference base from reference sequence
        ref_base = ref_sequence[pos]

        # Loops through each base and its count at the position
        for base, count in bases.items():
            frequency = count / total_coverage  # Calculates frequency of the base

            # Compares base with reference base to identify mutations
            if frequency != 1.0 and base != ref_base:
                mut_freq[pos] = (base, frequency * 100)  # Stores mutation information

    return mut_freq  # Returns dictionary containing mutation frequencies

def read_ref_seq(ref_file): # Reads reference sequence FASTA file to index bases
    with open(ref_file, 'r') as fasta_file: # Opens reference sequence FASTA file
        lines = fasta_file.readlines() # Reads lines in FASTA file
        ref_sequence = ''.join(line.strip() for line in lines[1:]) # Ignores first line, creates reference sequence string
    return ref_sequence


# Gets mold color information from data file
def read_data_file(file_path):
    color_dict = {} # Initializes dictionary to map names to colors   
    with open(file_path, 'r') as file: # Opens data file for reading
        lines = file.readlines() # Reads all lines from the file
        header = lines[0].strip().split('\t') # Gets header line and splits into fields using tab as delimiter
        for line in lines[1:]: # Loops through lines starting from second line (skip header)
            name, color, barcode = line.strip().split('\t') # Splits current line into fields using tab as delimiter
            color_dict[name] = color  # Maps name to corresponding color in dictionary
    return color_dict

def main():
    data_file = 'harrington_clinical_data.txt' # Specifies path to data file for color-name mapping
    bam_dir = 'BAM' # Specifies path to BAM directory
    ref_file = 'dgorgon_reference.fa' # Specifies path to reference sequence
    ref_sequence = read_ref_seq(ref_file) # Reads reference sequence FASTA into string
        
    color_name_mapping = read_data_file(data_file) # Reads data file to map sample names to colors
    
    with open('report.txt', 'w') as report_file: # Opens report.txt to create report file
        bam_files = [filename for filename in os.listdir(bam_dir) if filename.endswith("_trimmed.sorted.bam")]# Lists of BAM files that end with "_trimmed.sorted.bam"
        sorted_bam_files = sorted(bam_files) # Sorts list of BAM files alphabetically
        
        for filename in sorted_bam_files: # Loops through each sorted BAM file           
            sample_name = filename.replace("_trimmed.sorted.bam", "") # Extracts sample name from filename            
            color = color_name_mapping.get(sample_name, "Unknown") # Gets color associated with sample name, or set it as "Unknown" if not found      
            bam_file_path = os.path.join(bam_dir, filename) # Constructs full path to BAM file  
            ntdict, total_reads = pileup(bam_file_path) # Calculates nucleotide frequencies 
            mut_freq = calc_mut_freq(ntdict, total_reads,ref_sequence) # Calculates mutation frequencies 
            report_file.write(f"Sample {sample_name} had a {color} mold.\n") # Writes sample information and calculations to report file
            for pos, (base, frequency) in mut_freq.items(): 
                report_file.write(f"At position {pos+1}, {frequency:.2f}% of the reads had the mutation {base}.\n")
            
            report_file.write("\n")  

    print("Mutations and their frequencies saved to report.txt.")

if __name__ == "__main__":
    main()


