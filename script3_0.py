import sys

"""
function_placeholder
Author: Delshad Vegter, Bram Koobs 
date: 01/11/2022
version: 3
"""
def extract_microarray_content(gene_probes_dictionary):
    """
    Extract the highest value probes from the microarray for use
    :param gene_probes_dictionary: dict of lists
    :return: list of lists
    """

    # predefined variables
    average = 0
    microarray_data = {}
    index = 0
    highest_value_probe_gene = []
    highvalue_probes = []

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
        best_probes_list = probe_numbers[best_probe_index]

        # make a list with all the highest value probes
        highest_value_probe_gene.append(best_probes_list)
        best_probes_list = ""


    for lines in file_content:
        for probes in highest_value_probe_gene:
            if lines.startswith(probes):
                print(probes)
                highvalue_probes.append(lines)

    return highvalue_probes


def extract_sample_annot(identifier):
    """
    :returns 2 lists, each containing the corresponding line number from the requested structure_acronyms.
    """

    # predefined variables
    lines_identifier1 = []
    lines_identifier2 = []


    with open("./data/SampleAnnot.csv", "r") as file:
        sample_annot = list(file.readlines())

    for id in identifier:
        count = -1
        for line in sample_annot:
            count += 1
            if id in line:
                if id == identifier[0]:
                    lines_identifier1.append(count)
                else:
                    lines_identifier2.append(count)

    line_identifiers = [lines_identifier1, lines_identifier2]
    return line_identifiers

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

def extract_above_cutoff(line_identifiers, high_value_probes, cutoff_value):

    # predefined variables
    identifier_values1 = []
    identifier_values2 = []

    for lines in high_value_probes:
        lines = lines.split(",")
        probe_id = lines.pop(0)

        for index in line_identifiers[0]:
            if float(lines[int(index)]) > float(cutoff_value):
                identifier_values1.append(probe_id)

        for index in line_identifiers[1]:
            if float(lines[int(index)]) > float(cutoff_value):
                identifier_values2.append(probe_id)

    return identifier_values1, identifier_values2

def gene_name_finder(identifier_values1, identifier_values2):

    identifier_1_gene = []
    identifier_2_gene = []

    with open("./data/Probes.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    for lines in file_content:
        for values in identifier_values1:
            if lines.startswith(values):
                splitline = lines.split(",")
                identifier_1_gene.append(splitline[3])

        for values2 in identifier_values2:
            if lines.startswith(values2):
                splitline = lines.split(",")
                identifier_2_gene.append(splitline[3])

    return identifier_1_gene, identifier_2_gene

def main():
    """
    Main, where all functions are run in order
    :return: /
    """

    arguments = sys.argv

    if len(arguments) < 5:
        raise Exception(
            "No input filename given, please use commandline arguments. "
            "(python3 script2_3.py identifier identifier2)")
    else:
        # get the argument and set the variable
        identifier = [arguments[2], arguments[3]]
        cutoff_value = arguments[4]

    lines_identifiers = extract_sample_annot(identifier)
    gene_probes_dictionary = extract_genes()
    highvalue_probes = extract_microarray_content(gene_probes_dictionary)
    identifier_values1, identifier_values2 = extract_above_cutoff(lines_identifiers, highvalue_probes, cutoff_value)
    identifier_1_gene, identifier_2_gene = gene_name_finder(identifier_values1, identifier_values2)

    print(
        f'{len(identifier_1_gene)} unique genes for {arguments[2]}: {identifier_1_gene} \n \n{len(identifier_2_gene)} unique genes for {arguments[3]}: {identifier_2_gene}')

# to protect against problems when imported
if __name__ == "__main__":
    main()