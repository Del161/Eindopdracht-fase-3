import sys

"""
function_placeholder
Author: Delshad Vegter 
date: 01/11/2022
version: 1.1
"""
def extract_microarray_content(gene_probes_dictionary):
    """
    Opening of the given file, extract microarray data to dictionary.
    :return: list with content of the file
    """

    # predefined variables
    average = 0
    microarray_data = {}
    index = 0
    best_probes_list = []



    with open("./data/MicroarrayExpression.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    # take each line, turn the id into a key and the content into the key content.
    for line in file_content:
        line = line.strip().split(",")
        probe_id = line.pop(0)
        probe_content = line[1:]
        microarray_data[probe_id] = probe_content
        index += 1

    for genes in gene_probes_dictionary:
        probe_numbers = gene_probes_dictionary[genes]
        for probes in probe_numbers:
            average = 0
            rowvalues = []

            # split the row at the commas to get the seperate values
            microarray_row = microarray_data[probes]

            # combine all the values in a probe
            for values in microarray_row:
                average += float(values)

            # get the index of the highest value probe
            average = average / len(microarray_row)
            rowvalues.append(average)
            highest_value = max(rowvalues)
            best_probe_index = rowvalues.index(highest_value)
        best_probes_list.append(probe_numbers[best_probe_index])

        print(genes, best_probes_list)
        best_probes_list = []





    return

def extract_sample_annot(identifier):


    with open("./data/SampleAnnot.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    index = 0
    probes_data = {}
    locations = []

    # take each line, turn the id into a key and the content into the key content.
    for lines in file_content[1:]:
        lines = lines.strip().split(",")
        probe_ID = lines.pop(0)
        probe_content = "".join(lines)
        probes_data[probe_ID] = probe_content

    for lines in file_content:
        # find which line contains the identifier
        if identifier[0] in lines:
            locations.append(index)
        index += 1

    return probes_data

def extract_genes():

    # predefined vars
    gene_id = ""
    probe_id = ""
    gene_probeid_list = ""
    genes_list =[]
    gene_probes_dictionary = {}

    with open("./data/Probes.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())
    for lines in file_content[1:]:
        lines = lines.split(",")
        gene_id = lines.pop(2)
        probe_id = lines.pop(0)

        if gene_id in gene_probes_dictionary:
            gene_probes_dictionary[gene_id].append(probe_id)
        else:
            gene_probes_dictionary[gene_id] = list([probe_id])

    return gene_probes_dictionary


def main():
    """
    Main, where all functions are ran in order
    :return: /
    """

    arguments = sys.argv

    if len(arguments) < 3:
        raise Exception(
            "No input filename given, please use commandline arguments. "
            "(python3 script2_3.py identifier identifier2)")
    else:
        # get the argument and set the variable
        identifier = [arguments[1], arguments[2]]

    probes_data = extract_sample_annot(identifier)
    gene_probes_dictionary = extract_genes()
    extract_microarray_content(gene_probes_dictionary)


# to protect against problems when imported
if __name__ == "__main__":
    main()