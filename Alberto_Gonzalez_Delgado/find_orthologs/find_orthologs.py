#!/usr/bin/env python3

'''
Script for finding orthologs using reciprocal BLAST hits.

Sript usage is as follows:
    ./find_orthologs.py -i1 <Input file 1> -i2 <Input file 2> -o <Output file name> –t <Sequence type – n/p> -c <Number of target sequences> -s <Number of high-scoring segment pairs> -org1 <Organism code 1> -org2 <Organism code 2> -keep <Boolean: Keep the temporary files> 

where "n" specifies a nucleotide sequence and "p" specifies a protein sequence.
'''
import argparse
import subprocess

def get_reciprocal_hits(file_one,file_two, output_file,input_sequence_type,target_seqs,max_hsps,org_code1,org_code2,keep):

    a = []                  #List a stores the first index headers (db values)
    a1 = []                 #List a1 stores the second index headers (query headers)
    b = []                  #List b stores the first index headers (db values)
    b1 = []
    output_list_1 = []
    output_list_2 = []      #Inversing the list in order to find common things bw output_list_1 and output_list_2
    output_list = []        #Contains only reciprocal blast hits
#
    if input_sequence_type == "n":  #For nucleotide(n)

        file_1_db = "makeblastdb -in " + str(file_one) + " -dbtype nucl -out tmp/db1"   #Creating db for input1
        file_1_db_subprocess = subprocess.check_output(file_1_db.split())

        result_1 = "blastn -db tmp/db1 -query " + str(file_two) + " -max_target_seqs " +str(target_seqs) +" -max_hsps " +str(max_hsps) +" -outfmt 6 -out tmp/output_file_1" #Output for query = input2 and db = input1
        result_1_subprocess = subprocess.check_output(result_1.split())

        file_2_db = "makeblastdb -in " + str(file_two) + " -dbtype nucl -out tmp/db2"   #Creating db for input2
        file_2_db_subprocess = subprocess.check_output(file_2_db.split())

        result_2 = "blastn -db tmp/db2 -query " + str(file_one) + " -max_target_seqs " +str(target_seqs) +" -max_hsps " +str(max_hsps) +" -outfmt 6 -out tmp/output_file_2" #Output for query = input1 and db = input2
        result_2_subprocess = subprocess.check_output(result_2.split())


        with open("tmp/output_file_1","r") as fh1:
            for line in fh1:
                if line.startswith("lcl"):
                    a.append(line.split()[0])           #List a stores the first index headers (db values)
                    a1.append(line.split()[1])          #List a1 stores the second index headers (query headers)

        for i,j in zip(a,a1):
            output_list_1.append(i+"\t"+j+"\n")

        with open("tmp/output_file_2","r") as fh2:
            for line in fh2:
                if line.startswith("lcl"):
                    b.append(line.split()[0])           #List b stores the first index headers (db values)
                    b1.append(line.split()[1])          #List b1 stores the second index headers (query headers)

        for i,j in zip(b,b1):
            output_list_2.append(j+"\t"+i+"\n")         #Inversing the list in order to find common things bw output_list_1 and output_list_2

        for x in output_list_1:                         #Checking common hits beween output_list_1 and output_list_2
            if x in output_list_2:
                output_list.append(x)                   #Contains only reciprocal blast hits

    elif input_sequence_type == "p":  #For protein(p)

        file_1_db = "makeblastdb -in " + str(file_one) + " -dbtype prot -out tmp/db1"   #Creating db for input1
        file_1_db_subprocess = subprocess.check_output(file_1_db.split())
#max target and max_hsps 1

        result_1 = "blastp -db tmp/db1 -query " + str(file_two) + " -max_target_seqs " +str(target_seqs) +" -max_hsps " +str(max_hsps) +" -outfmt 6 -out tmp/output_file_1" #Output for query = input2 and db = input1
        result_1_subprocess = subprocess.check_output(result_1.split())

        file_2_db = "makeblastdb -in " + str(file_two) + " -dbtype prot -out tmp/db2"   #Creating db for input2
        file_2_db_subprocess = subprocess.check_output(file_2_db.split())

        result_2 = "blastp -db tmp/db2 -query " + str(file_one) + " -max_target_seqs " +str(target_seqs) +" -max_hsps " +str(max_hsps) +" -outfmt 6 -out tmp/output_file_2" #Output for query = input1 and db = input2
        result_2_subprocess = subprocess.check_output(result_2.split())

        with open("tmp/output_file_1","r") as fh1:
            for line in fh1:
                if line.startswith(f'{org_code1}'):
                    a.append(line.split()[0])
                    a1.append(line.split()[1])

        for i,j in zip(a,a1):
            output_list_1.append(i+"\t"+j+"\n")

        with open("tmp/output_file_2","r") as fh2:
            for line in fh2:
                if line.startswith(f'{org_code2}'):
                    b.append(line.split()[0])
                    b1.append(line.split()[1])

        for i,j in zip(b,b1):
            output_list_2.append(j+"\t"+i+"\n")
        for x in output_list_1:         #Checking common hits beween output_list_1 and output_list_2
            if x in output_list_2:
                output_list.append(x)   #Contains only reciprocal blast hits

    return(output_list)

def main():
    # Argparse code

    parser = argparse.ArgumentParser()
    parser.add_argument("-i1", help = "Takes in input sequence 1.")
    parser.add_argument("-i2", help = "Takes in input sequence 2.")
    parser.add_argument("-o", help = "Your output file.")
    parser.add_argument("-t", help = "n for nucleotide sequence, p for protein sequence.")
    parser.add_argument("-c", help = "target secuences for comparing (recommended at least 5)")
    parser.add_argument("-s", help = "max number of High-scoring Segment pair.")
    parser.add_argument("-org1", help = "Organism code for input 1")
    parser.add_argument("-org2", help = "Organism code for input 2")
    parser.add_argument("-keep", help = "Keep the temporary files", action='store_true')
    args = parser.parse_args()


    file_one = args.i1             
    file_two = args.i2
    output_file = args.o
    input_sequence_type = args.t
    target_seqs = args.c
    max_hsps = args.s
    org_code1 = args.org1
    org_code2 = args.org2
    keep=args.keep

    '''
    output_list is a list of reciprocal BLAST hits. Each element is a tab
    separated pair of gene names. The list is written to the output file.
    '''

    make_dir = "mkdir tmp"
    make_file_1='touch tmp/output_file_1'
    make_file_2='touch tmp/output_file_2'

    subprocess.call(make_dir.split())
    subprocess.call(make_file_1.split())
    subprocess.call(make_file_2.split())

    output_list = get_reciprocal_hits(file_one, file_two, output_file, input_sequence_type, target_seqs, max_hsps, org_code1, org_code2,keep)
    print(output_list)
    with open(output_file, 'w') as output_fh:
        for ortholog_pair in output_list:
            output_fh.write(ortholog_pair)
    if keep:
        copy_file_1='cp tmp/output_file_1 ./blast_output_1.tsv'
        copy_file_2='cp tmp/output_file_2 ./blast_output_2.tsv'
        subprocess.call(copy_file_1.split())
        subprocess.call(copy_file_2.split())

    remove_dir = "rm -rf tmp"
    subprocess.call(remove_dir.split())

if __name__ == "__main__":
    main()
