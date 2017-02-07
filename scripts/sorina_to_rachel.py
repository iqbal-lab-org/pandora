import sys
from Bio import SeqIO

input_file = sys.argv[1]
output_file = "".join(input_file.split(".")[:-1]) + ".rachel.fa"
print input_file
print output_file

fasta_sequences = SeqIO.parse(open(input_file,'r'),'fasta')
with open(output_file, 'w') as out_file:
    for fasta in fasta_sequences:
        name, description, sequence = fasta.id, fasta.description, str(fasta.seq)
	new_sequence = ""
	num = ""
	for letter in sequence:
	    if letter.isdigit():
		num += letter
	    elif num!="":
		new_sequence += " " + num + " " + letter
		num = ""
	    else:
		new_sequence += letter
        out_file.write(">%s\t%s\n%s\n" %(name, description, new_sequence))
