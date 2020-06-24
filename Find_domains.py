#!/usr/bin/env python

import re

from Bio.ExPASy import Prosite,Prodoc
from Bio import SeqIO

def Parsear_prosite():

    """Function to parse the data base selecting the elements
    necessary for the create file"""

    with open("prosite.dat", "r") as data_base,\
            open("prosite_parser.tsv", "w") as result_file:
        result_file.write("Name"+"\t"+"Accesion"+"\t"+"Desctiption"+"\t"+"Pattern"+"\n")
        records = Prosite.parse(data_base)

        for record in records:
            name = record.name
            accesion = record.accession
            description = record.description
            patt = record.pattern
            result = str(name+"\t"+accesion+"\t"+description+"\t"+patt+"\n")
            result_file.write(result)

        records.close()
        result_file.close()

def Convert_regex (file_pros='prosite_parser.tsv'):

    """Funtion for being able to use the regular expresion
    and serach with them"""

    final_patron = file_pros.replace(".", "")
    final_patron = final_patron.replace("{", "[^").replace("}", "]")
    final_patron = final_patron.replace("(", "{").replace(")", "}")
    final_patron = final_patron.replace("<", "^").replace(">", "$")
    final_patron = final_patron.replace("x", ".")
    final_patron = final_patron.replace("-", "")

    return (final_patron)

def find_domains (file_query, folder_result, file_pros="prosite_parser.tsv"):

    """Function to find the protein domain of our file in the data base.
    We select to write name,a ccesion, description and pattern. THe process will
    be done por each query and save in the correct folder"""

    for record in SeqIO.parse(file_query, "fasta"):

        with open(record.id+'_muscle.faa', "r") as result:
            with open(file_pros, "r") as file_prosite,\
                    open(record.id+"_domains.txt", "w") as finds:
                division_result = result.read().split("\n")
                division_prosite = file_prosite.read().split("\n")

                for i in range(len(division_result) // 2):
                    sequence = division_result[2*i + 1]
                    name = division_result[2*i]
                    seq = str("Pattern in:"+name[1:]+"\n\n")
                    finds.write(seq)

                    for j in division_prosite[1:]:
                        #Not to read the firt line, header
                        if j == '':
                            pass
                        else:
                            j = j.split("\t")

                            #Read pattern with the regular expresion

                            find_patron=Convert_regex(j[3])
                            if j[3] != "":
                                if re.search(find_patron, sequence):
                                    res = str("\tName"+j[0]+"\n\tAccesion:"+j[1]+""
                                            "\n\tDescription:"+j[2]+"\n\tPattern:"+j[3])
                                    finds.write(res)

            file_prosite.close()

    print("The domains that were found are save in the folder: %s" %(folder_result))

