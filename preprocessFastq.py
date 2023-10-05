import os
import argparse
import gzip

class ParseFastQ(object):
    
    def __init__(self,filePath,headerSymbols=['@','+']):
        
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
        
# Reads clinical data file using tab as delimiter; saves data to dictionary
def read_file(file_path):
    clin_data_dict = {} # Initializes dictionary to map names to barcodes 
    with open(file_path, 'r') as file: # Opens data file for reading
        lines = file.readlines() # Reads all lines from the file
        header = lines[0].strip().split('\t') # Gets header line and splits into fields using tab as delimiter
        for line in lines[1:]: # Loops through lines starting from second line (skip header)
            name, _, barcode = line.strip().split('\t')  # Splits current line into fields using tab as delimiter
            clin_data_dict[barcode] = name # Maps name to corresponding barcode in dictionary
    return clin_data_dict
    
# Deletes barcode from read sequence
def del_barcode(sequence, barcode):
    if sequence.startswith(barcode): # Checks if read sequence starts with barcode
        return sequence[len(barcode):] # If sequence starts with barcode, returns part of sequence after barcode
    return sequence # If sequence doesn't start with barcode, returns sequence as is

# Trims end of reads if has degraded quality score of at least 2 consecutively
def trim_read_ends(sequence, quality_score):
    consecutive_count = 0 # Counter tracks consecutive degraded quality scores
    for i, q in enumerate(quality_score):  # Checks if quality score is degraded base
        if q == 'D' or q == 'F':
            consecutive_count += 1 # Increases count of consecutive degraded scores
            if consecutive_count >= 2:  # If at least two consecutive D or F quality scores
                return sequence[:i-1], quality_score[:i-1]  # Removes corresponding portion of sequence and quality scores
        else:
            consecutive_count = 0  # Resets counter if good quality score comes next
    return sequence, quality_score # If consecutive degraded quality scores not found, returns original sequence and quality scores

def process_fastq(fastq_file, barcode_name_mapping):
  
    # Loop through each fastq_obj read from the ParseFastQ instance
    for fastq_object in fastq_file:
        seq_name = fastq_object[0]
        sequence = fastq_object[1]
        separator = fastq_object[2]
        quality_score = fastq_object[3]

        # Find the matching barcode and corresponding name
        matched_name = None
        for barcode, name in barcode_name_mapping.items():
            if sequence.startswith(barcode):
                matched_name = name
                break
        
        if matched_name is None:
            continue  # Skip sequence if no matching barcode was found
        
        # Trims sequence by deleting barcode and trims sequence ends by quality score
        trimmed_sequence = del_barcode(sequence, barcode)
        trimmed_sequence, trimmed_quality_score = trim_read_ends(trimmed_sequence, quality_score)

        # Writes fastq sequence to corresponding matched name file
        trimmed_file_path = f'fastqs/{matched_name}_trimmed.fastq'
        with open(trimmed_file_path, 'a') as trimmed_file:
            trimmed_file.write(seq_name + '\n')
            trimmed_file.write(trimmed_sequence + '\n')
            trimmed_file.write(separator + '\n')
            trimmed_file.write(trimmed_quality_score + '\n')

def main():
    if not os.path.exists('fastqs'): # Creates fastqs director if none exists
        os.makedirs('fastqs')

    parser = argparse.ArgumentParser()  # Parse command line arguments to get the path of the input FASTQ file
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    args = parser.parse_args()

    file_path = 'harrington_clinical_data.txt' # Path to file containing data for barcode-to-name mapping
    barcode_name_mapping = read_file(file_path) # Maps barcodes to names 

    fastqfile = ParseFastQ(args.fastq) # Initializes ParseFastQ instance to read input pooled FASTQ sequence file

    process_fastq(fastqfile, barcode_name_mapping) # Processes each FASTQ read sequence 

if __name__ == "__main__":
    main()

