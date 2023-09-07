from Bio import SeqIO
import os
import sys

data = open(sys.argv[3], "w+")
with open(sys.argv[1]) as f:
    for id in f:
        print(id.rstrip())
        seq = SeqIO.parse(sys.argv[2], "fasta")
        for line in seq:
            if line.id == id.rstrip():

                data.write(f">{str(id.rstrip())}" +
                           "\n" + str(line.seq) + "\n")
