#!/usr/bin/env python
#-*- coding: UTF-8 -*_

import os
import pandas as pd
import sys

from Bio import SeqIO

from subprocess import Popen, PIPE


def Parser(folder):

    """ Funtion to determined if the files have the format genbank
    and to rewrite the file in fasta format """

    directory = os.getcwd()

    try:
        if os.path.isfile("file_parseado") == True:
            os.remove("file_parseado")
        else:
            pass
        output_handle = open("file_parseado", "a")
        os.chdir(folder)
        os.listdir()
        for file in os.listdir():
            with open(file, "r") as input_handle:
                for seq_record in SeqIO.parse(input_handle, "genbank"):
                    for seq_feature in seq_record.features:
                        try:
                            if seq_feature.type == 'CDS':

                                #Select the outputs for the file

                                output_handle.write("> %s@%s\n%s\n" % (
                                    seq_feature.qualifiers['locus_tag'][0],
                                    seq_record.name,
                                    seq_feature.qualifiers['translation'][0]))
                        except:
                            pass
            input_handle.close()
        output_handle.close()
        print("Done!")
        os.chdir(directory)
    except:
        print("Incorrect format")
        sys.exit()


def BlastP(file_query, file_parseado):

    """Funtion for doing blastp where is needed two differents files.
    Files must be fasta. The result will be save as {id}_result.tsv
    in the folder choose by the user"""

    with open(file_query, "r") as query:
        for record in SeqIO.parse(query, "fasta"):
            seq = str(">%s\n%s" % (record.id, record.seq))
            temp = open(record.id + "temp.faa", "w")
            temp.write(seq)
            temp.close()
            # Call the function blast
            with open(record.id + "_result.tsv", "w") as res_blast:
                process = Popen(['blastp', '-query', record.id + "temp.faa", '-subject',
                                 file_parseado, '-evalue', '0.00001', '-outfmt',
                                 "6 qseqid qcovs pident evalue sseqid sseq"], stdout=PIPE, stderr=PIPE)
                header = str("query_ID\tCoverage\tIdentity\tEvalue\tSubject_ID\tSubject_seq\n")
                result = process.stdout.read().decode("utf-8")
                res_blast.write(header)
                res_blast.write(result)
                res_blast.close()
                os.remove(record.id + "temp.faa")

        print("Blastp is done")


def Values(file_query):

    """Funtion to filter the blastp result. First we save the values of
     coverage and identity in variables. In case they are not introduced
     the variable will be assign with exactly values. The the file will
     be filter. Outputs will be save in the folder created and in format
     tsv"""

    # Control the value of coverage

    try:
        if int(sys.argv[3]) >= 100:
            print("Error: coverage has to be in between 0 and 100. Try again.")
            sys.exit()
        elif int(sys.argv[3]) >= 0 or int(sys.argv[3]) <= 100:
            print("This coverage is OK.")
            value_cover = sys.argv[3]
            pass
    except:
        value_cover = 50

    # Control de value of identity

    try:
        if int(sys.argv[4]) >= 100:
            print("Error: identity has to be in between 0 and 100. Try again.")
            sys.exit()
        elif int(sys.argv[4]) >= 0 or int(sys.argv[4]) <= 100:
            print("This identity is OK.")
            value_identity = sys.argv[4]
            pass
    except:
        value_identity = 30

    # Filtration

    for record in SeqIO.parse(file_query, "fasta"):
        with open(record.id + "_result.tsv") as tsvfile, \
                open(record.id + "_filter.tsv", "w") as tsv_filter:
            tsvreader = pd.read_csv(tsvfile, delimiter='\t')
            trying = tsvreader.loc[(tsvreader['Identity'] >= int(value_identity))
                                   & (tsvreader['Coverage'] >= int(value_cover)), :]
            trying.to_csv(tsv_filter, sep='\t')
            tsvfile.close()
            tsv_filter.close()
    print("The filtration is done correctly")