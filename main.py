#!/usr/bin/env python
#-*- coding: UTF-8 -*_

import os
import shutil
import sys

from Bio import SeqIO

#Importar the function in other script
import Find_domains
import funtion_blast
import make_tres

#Informatatio for teh program help
if len(sys.argv)<2:
    print("""\n --------------Welcome to the program help-------------
    
    The aim of this program is doing a blastp, a muscle aligment and 
    find the domains of the protein sequence. 
    For doing the blast will pass from the genbak format to a fasta format
    for doing the blastp. Having one file for each query_ID. Then the
    result will be filtered with the value introduced by the user. If the 
    user dont introduced any value, the program has definied the values. 
    
    The muscle alignment is done with the querys and their subjects. With
    this output will can create the phylogenetic tree. 
    
    The domains where found using the file 'prosite.dat' which should
    be located in the folder. This output wil be used to find the pattern.
    All the results will be copy in a folder with the name that the user 
    has choose. 
    
    To run the program you hace to introduced in this order:
        1. The name of folder where the the files of the genbanks
            are save.
        2. Query file
        3. Value of coverage (not compolsary)
        4. Value of identity (not compolsary
    
    --------------Thank ok for using the help--------------
    """)
    sys.exit()
else:
    pass
#Folder to save the results
print("Introduced name for the folder to save the results")
usr_ip=input()
try:
    folder_result="Resultado_" + str(usr_ip)
    os.mkdir(folder_result)
except OSError:
    print("Error creating the folder. Please try again")
    exit()
else:
    print("Succesful creating the directory")

#Assigment of value

folder=sys.argv[1]
file_query=sys.argv[2]

#Call to the function Parser

funtion_blast.Parser(folder)
shutil.copy('file_parseado', folder_result)

#Call to the function blast

funtion_blast.BlastP(file_query, 'file_parseado')

#Exception when the values are not correct

funtion_blast.Values(file_query)

#Call the function for alignment

make_tres.Pre_Alignment(file_query)
make_tres.Alignment(file_query, 'folder_result')
make_tres.Tree(file_query)

#Call to the function to find domains

Find_domains.Parsear_prosite()
Find_domains.find_domains(file_query,folder_result, file_pros='prosite_parser.tsv')

#Move the file to the correct folder

for record in SeqIO.parse(file_query, "fasta"):
    shutil.copy(record.id + "_result.tsv", folder_result)
    shutil.copy(record.id + "_filter.tsv", folder_result)
    shutil.copy(record.id + "_muscle.faa", folder_result)
    shutil.copy(record.id + "_alignment.faa", folder_result)
    shutil.copy(record.id + "_tree.nw", folder_result)
    shutil.copy(record.id+"_domains.txt", folder_result)
