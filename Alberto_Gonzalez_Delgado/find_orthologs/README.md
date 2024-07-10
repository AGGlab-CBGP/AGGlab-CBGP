# Find Orthologs
Reciprocal orthologs lists (-c number of reciprocal).

## Requirements:
BLAST+ for running makeblastdb and blast.

## Required input files
1. FASTA file for specie 1.
2. FASTA file for specie 2.

## Usage:
    ./find_orthologs.py -i1 <Input file 1> -i2 <Input file 2> -o <Output file name> –t <Sequence type – n/p> -c <Number of target sequences> -s <Number of high-scoring segment pairs> -org1 <Organism code 1> -org2 <Organism code 2> -keep <Optional: Keep the temporary files if given> 
Organism code refers to start of gene_id names in FASTA files.
## Output
1. A file comprising the reciprocal BLAST hits of the two species in output format 6.
2. A file comprising the BLAST results for Organism 1.
3. A file comprising the BLAST results for Organism 2.



