#!/usr/bin/env python
# -- coding: utf-8 --

from Bio import SeqIO

from subprocess import Popen, PIPE

def Pre_Alignment (file_query, blast_filter="blastp_filter.tsv"):

    """Funtion for preparing the files for the aligment with
    muscle. THe files will be in format fasta withthe subject
    ID and the sequence. File will be saved in the folder result"""

    for record in SeqIO.parse(file_query, "fasta"):
        seq=str(">%s\n%s\n" %(record.id, record.seq))
        test = open(record.id+"_filter.tsv", "r")
        out_muscle = open(record.id+"_muscle.faa", "w")
        out_muscle.write(seq)

        with open(record.id+"_muscle.faa", "a+") as out_muscle:
            for row in test.readlines():
                fields=row.rstrip().split("\t")
                if fields[6] !='Subject_seq':
                    Subject_ID = fields[5]
                    Subject_seq = fields[6]
                    out_muscle.write(">%s\n%s\n" % (Subject_ID, Subject_seq))

        out_muscle.close()

def Alignment (file_query, folder_result):

    """Function to alignment using muscle. File for each query in
     fasta format will be used. The output will be a file for
     each query with the alignment done. Also will be save in the
     folder result"""

    for record in SeqIO.parse(file_query, "fasta"):
        with open(record.id+"_muscle.faa", "r") as file,\
                open(record.id+"_alignment.faa", "a") as alignment:
            process = Popen(['muscle', '-in', record.id + "_muscle.faa", "-out",
                             record.id + "_alignment.faa"], stdout=PIPE, stderr=PIPE)
            result = process.stdout.read().decode("utf-8")
            process.stdout.close()
            alignment.write(result)
        alignment.close()
        file.close()

    print ("Alignment is done")

def Tree (file_query, input_alig='record.id+"_alignment.faa"'):

    """Function for doing a neighborjoining tree. There will be a diferent file
    for each query"""

    for record in SeqIO.parse(file_query, "fasta"):
        with open(record.id + "_alignment.faa", "r") as input_alig,\
                open(record.id + "_tree.nw", "w") as output_tree:
            tree = Popen(['muscle', '-maketree', '-in', record.id + "_alignment.faa",
                          '-out', record.id + "_tree.nw",'-cluster', 'neighborjoining'],
                         stderr=PIPE)
        input_alig.close()
        output_tree.close()

    print("Neighborjoining tree was made correctly")