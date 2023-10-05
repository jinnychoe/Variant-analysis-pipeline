import preprocessFastq  # Import module for preprocessing FASTQ files
import alignFastq       # Import module for aligning FASTQ files
import samToBam         # Import module for converting SAM to BAM files
import findVariant     # Import module for finding variants in aligned sequences

def main():
    preprocessFastq.main()  # Execute main function in preprocessFastq module
    alignFastq.main()       # Execute main function in alignFastq module
    samToBam.main()         # Execute main function in samToBam module
    findVariant.main()     # Execute main function in findVariants module

if __name__ == "__main__":
    main()  
