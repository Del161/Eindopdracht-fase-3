import sys
import os.path # to check if a improved file already exists
"""
function_placeholder
Author: Delshad Vegter, Bram Koobs 
date: 01/11/2022
version: 4

To use enter in the commandline: python3 script3_0.py (wanted_symbol) identifier1 identifier2 cutoff_value
example: python3 script3_0.py gene_symbol LHM PHA 17
possible wanted symbols: gene_symbol, gene_id, gene_name, entrez_id, chromosome
"""


def extract_microarray_content(gene_probes_dictionary, identifiers):
    """
    Extract the highest value probes from the micro-array for use
    :param gene_probes_dictionary: dict of lists
    :return: list of lists
    """

    # predefined variables
    average = 0
    microarray_data = {}
    highest_value_probe_gene = []
    high_value_probes = []
    rowvalues = []

    with open("./data/MicroarrayExpression.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    # take each line, turn the id into a key and the content into the key content.
    for line in file_content:
        line = line.strip().split(",")
        probe_id = line.pop(0)
        probe_content = line
        microarray_data[probe_id] = probe_content

    print("created dictionary...")

    for genes in gene_probes_dictionary:
        probe_numbers = gene_probes_dictionary[genes]
        for probes in probe_numbers:
            average = 0

            microarray_row = microarray_data[probes]

            # combine all the values in a probe
            for values in microarray_row:
                average += float(values)

            # get the index of the highest value probe
            average = average / len(microarray_row)
            rowvalues.append(average)

        # get the index of the probe with the highest value, and append it to a list
        highest_value = max(rowvalues)
        best_probe_index = rowvalues.index(highest_value)
        best_probes_list = probe_numbers[best_probe_index]
        rowvalues = []
        # make a list with all the highest value probes
        highest_value_probe_gene.append(best_probes_list)
        best_probes_list = ""

    print("got highest values...")

    for lines in file_content:
        if lines.startswith(tuple(highest_value_probe_gene)):
            high_value_probes.append(lines)
    print("made high value list...")

    with open("High_Value_Probes.csv", "w") as writefile:
        for lines in file_content:
            if lines.startswith(tuple(highest_value_probe_gene)):
                writefile.write(lines)

    return high_value_probes

def use_premade_probes():

    with open("High_Value_Probes.csv", "r") as file:
        # set file lines to variable
        high_value_probes = list(file.readlines())

    return  high_value_probes


def extract_sample_annot(identifier):
    """
    :param list of strings
    :returns 2 lists, each containing the corresponding line number from the requested structure_acronyms.
    """

    # predefined variables
    lines_identifier1 = []
    lines_identifier2 = []

    with open("./data/SampleAnnot.csv", "r") as file:
        sample_annot = list(file.readlines())

    for ID in identifier:
        count = 0
        for line in sample_annot[1:]:
            count += 1
            if ID in line:
                if ID == identifier[0]:
                    lines_identifier1.append(count)
                else:
                    lines_identifier2.append(count)

    line_identifiers = [lines_identifier1, lines_identifier2]
    return line_identifiers


def extract_genes():
    """
    Extract the genes from the file, and put them in a dictionary with the probe id as key
    :return: dict of strings with lists.
    """

    # predefined vars
    gene_probes_dictionary = {}

    with open("./data/Probes.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    for lines in file_content[1:]:
        # split the lines, take the gene id and probe id
        lines = lines.split(",")
        gene_id = lines[3]
        probe_id = lines[0]

        # create a dictionary of gene id keys with probe id content
        if gene_id in gene_probes_dictionary:
            gene_probes_dictionary[gene_id].append(probe_id)
        else:
            gene_probes_dictionary[gene_id] = list([probe_id])

    return gene_probes_dictionary


def extract_above_cutoff(line_identifiers, high_value_probes, cutoff_value):
    """
    extract the values at the correct index, if they are
    :param line_identifiers:
    :param high_value_probes:
    :param cutoff_value:
    :return: 2 lists
    """

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


def gene_name_finder(identifier_values1, identifier_values2, requested_symbol):
    identifier_1_gene = []
    identifier_2_gene = []

    with open("./data/Probes.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    if requested_symbol == "gene_id":
        requested_index = 2
    elif requested_symbol == "gene_symbol":
        requested_index = 3
    elif requested_symbol == "gene_name":
        requested_index = 4
    elif requested_symbol == "entrez_id":
        requested_index = 5
    elif requested_symbol == "chromosome":
        requested_index = 6

    # check if the value at the indicated index is over the cutoff value, if it is append the probe id
    for lines in file_content:
        for values in identifier_values1:
            if lines.startswith(values):
                splitline = lines.split(",")
                identifier_1_gene.append(splitline[int(requested_index)].strip())

        # repeat for the second given identifier
        for values2 in identifier_values2:
            if lines.startswith(values2):
                splitline = lines.split(",")
                identifier_2_gene.append(splitline[int(requested_index)].strip())

    return set(identifier_1_gene), set(identifier_2_gene)


def unique_or_shared(identifier_1_gene, identifier_2_gene):
    """
    take the extracted genes, and create new lists with the shared genes, and unique genes.
    :param identifier_1_gene: list of genes for identifier 1
    :param identifier_2_gene: list of genes for identifier 2
    :return: list of shared genes, list of unique genes 1 and 2
    """

    # prepare lists for use
    identifier_1_gene = set(identifier_1_gene)
    identifier_2_gene = set(identifier_2_gene)
    identifier_1_gene = list(identifier_1_gene)
    identifier_2_gene = list(identifier_2_gene)
    common_genes = []
    unique_genes_1 = []
    unique_genes_2 = []

    # extract the genes that appear in both lists
    if len(identifier_1_gene) > len(identifier_2_gene):
        # check if the genes in the shorter list are in the longer list
        for genes in identifier_2_gene:
            if identifier_1_gene.count(genes) > 0:
                common_genes.append(genes)
    else:
        for genes in identifier_1_gene:
            if identifier_2_gene.count(genes) > 0:
                common_genes.append(genes)

    # extract the genes that are unique to each list
    for genes in identifier_1_gene:
        if genes not in common_genes:
            unique_genes_1.append(genes)
    for genes in identifier_2_gene:
        if genes not in common_genes:
            unique_genes_2.append(genes)

    return common_genes, unique_genes_1, unique_genes_2


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
        requested_symbol = arguments[1]

    print("getting identifiers...")
    lines_identifiers = extract_sample_annot(identifier)
    print("extracting genes...")
    gene_probes_dictionary = extract_genes()
    print("calculating highest value probes...")
    if os.path.isfile("High_Value_Probes.csv"):
        high_value_probes = use_premade_probes()
    else:
        high_value_probes = extract_microarray_content(gene_probes_dictionary, identifier)
    print("removing values below cutoff...")
    identifier_values1, identifier_values2 = extract_above_cutoff(lines_identifiers, high_value_probes, cutoff_value)
    print("extracting gene names...")
    identifier_1_gene, identifier_2_gene = gene_name_finder(identifier_values1, identifier_values2, requested_symbol)
    print("checking similarities between genes...")
    common_genes, unique_genes_1, unique_genes_2 = unique_or_shared(identifier_1_gene, identifier_2_gene)
    print("done!")

    print(
        f'{len(unique_genes_1)} unique genes for {arguments[2]}: {unique_genes_1} \n{len(unique_genes_2)} unique genes for {arguments[3]}: {unique_genes_2}')
    print(len(common_genes), "shared genes", common_genes)


# to protect against problems when imported
if __name__ == "__main__":
    main()
