
# script para parsear la base de datos prosite presentes en los archivos
# prosite.doc y prosite.dat utilizando el modulo Biopython

from Bio.ExPASy import Prosite,Prodoc

# con este script podeis parsear el archivo .dat
# handle = open("prosite.dat","r")
# records = Prosite.parse(handle)
# for record in records:
# 	print("name:"+record.name)
# 	print("accession:"+record.accession)
# 	print("description:"+record.description)
# 	print("pattern:"+record.pattern)

# con este script podemos parsear el archivo .doc
handle = open("prosite.doc")
records = Prodoc.parse(handle)
for record in records:
	print(record.accession)
	#print(record.prosite_refs)
	print(record.text)
	print(count)
	#print(record.references)
