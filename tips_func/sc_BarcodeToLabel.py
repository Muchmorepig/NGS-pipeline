# python 3
import click  # pip install click

# Dictionary
dict_sequence_label = {}
# Key = cell sequence ; value : cell label
# Opening file with the conversion between sequence and cell label
# Column 1 : Single-Cell Barcode; Column 2 : Label (CD4, B, CD8 etc...)
input_directory = click.prompt(
    '> Enter absolute path of input files directory')
output_directory = click.prompt('> Enter absolute path of output directory')

with open(input_directory + "Convert_UMI_Label.tsv", "r") as file_cell_label:
    next(file_cell_label)
    #useless first line
    for line in file_cell_label:
        clean_line = line.rstrip("\n").split("\t")
        sequence = clean_line[0].replace("-", ".")  #Store barcode
        cell_label = clean_line[1]  # Store Cell Label
        dict_sequence_label[sequence] = cell_label

print("dictionnary done !")

header = []
# Now i got the conversion between sequence and cell type
# i'll read the file with gene count per cell sequence
# and write a new file with the dictionnary above
with open(input_directory + "Gene_Count_Per_Cell.tsv", "r") as file_gene_count:
    first_line_cell_sequence = file_gene_count.readline()
    #The first line of export in R present a mistake
    #save it and correct after
    #next(file_gene_count) 
    #Delete first line
    with open(output_directory + "preSigMatrix.tsv", "w") as file_matrix:
        split_first_line_cell_sequence = first_line_cell_sequence.rstrip(
            "\n").split("\t")
        for element in split_first_line_cell_sequence:
            header.append(dict_sequence_label[element])
        file_matrix.write("GENES\t" + '\t'.join(header) + "\n")
        # Once I wrote all the cells label I can write the gene counts
        for line in file_gene_count:
            file_matrix.write(line)
# The script in finished and the file is now created.

print("Single Cell Raw Counts Annotated Matrix susccessfully created!")