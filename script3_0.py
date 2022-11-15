"""
Takes information from a multitude of files, and returns the information the user requested
Author: Delshad Vegter
date: 01/11/2022
version: 5
To use, enter into the commandline:
            python3 script3_0.py (wanted_symbol) identifier1 identifier2 cutoff_value
example: python3 script3_0.py gene_symbol LHM PHA 17
possible requested symbols: gene_symbol, gene_id, gene_name, entrez_id, chromosome
if the found gene does not have a chromosome or entrez id,-
-they get a replacement nan value to still show comparison.
"""
import sys
import os.path  # to check if an improved file already exists


def read_microarray_content():
    """
    If there already is a file with high value probes, this function will extract the data from it.
    :return: list of strings
    """
    with open("./data/MicroarrayExpression.csv", "r") as file:
        # set file lines to variable
        microarray_content = list(file.readlines())

    return microarray_content


def extract_microarray_content(gene_probes_dictionary, microarray_content):
    """
    Extract the highest value probes from the microarray for use
    :param microarray_content list of strings
    :param gene_probes_dictionary: dict of lists
    :return: list of lists
    """

    # predefined variables
    microarray_data = {}
    highest_value_probe_gene = []
    high_value_probes = []
    row_values = []

    # take each line, turn the id into a key and the content into the key content.
    for line in microarray_content:
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
            row_values.append(average)

        # get the index of the probe with the highest value, and append it to a list
        highest_value = max(row_values)
        best_probe_index = row_values.index(highest_value)
        row_values = []
        # make a list with all the highest value probes
        highest_value_probe_gene.append(probe_numbers[best_probe_index])

    print("got highest values...")

    # create a list with only the high value probes.
    for lines in microarray_content:
        if lines.startswith(tuple(highest_value_probe_gene)):
            high_value_probes.append(lines)
    print("made high value list...")

    return high_value_probes


def write_intermediate_file(high_value_probes):
    """
    Creates an intermediate file with high value probes
     if there isn't already one
    :param high_value_probes: list of strings
    :return: /
    """
    # create a file with all the high value probes, so it doesn't have to be repeated the next run.
    with open("High_Value_Probes.csv", "w") as writefile:
        for lines in high_value_probes:
            writefile.write(lines)


def use_premade_probes():
    """
    if a pre-made list with high value probes already exists, use that instead of creating one again.
    return: list of strings
    """

    with open("High_Value_Probes.csv", "r") as file:
        # set file lines to variable
        high_value_probes = list(file.readlines())

    return high_value_probes


def extract_sample_annot(identifier):
    """
    extracts important info from the sample annotation
    :param identifier list of strings
    :returns 2 lists, each containing the corresponding
    line number from the requested structure_acronyms.
    """

    # predefined variables
    lines_identifier1 = []
    lines_identifier2 = []

    with open("./data/SampleAnnot.csv", "r") as file:
        sample_annot = list(file.readlines())

    for ids in identifier:
        # count is 1 lower than it should be, because
        # the horizontal index starts at 0 not 1
        count = 0
        # this starts at 1 to skip the header
        for line in sample_annot[1:]:
            count += 1
            if ids in line:
                if ids == identifier[0]:
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


def read_probes_content():
    """
    reads the content of probes.csv
    :return: list of strings
    """
    with open("./data/Probes.csv", "r") as file:
        # set file lines to variable
        probes_file_content = list(file.readlines())
    return probes_file_content


def gene_name_finder(identifier_values1, identifier_values2, requested_symbol, probes_file_content):
    """
    Finds the requested genes, and makes them into a list
    :param probes_file_content: list of strings
    :param identifier_values1: list of strings
    :param identifier_values2: list of strings
    :param requested_symbol: string
    :return: 2 lists of strings
    """
    identifier_1_gene = []
    index_dict = {
        "gene_id": 2,
        "gene_symbol": 3,
        "gene_name": 4,
        "entrez_id": 5,
        "chromosome": 6
    }
    identifier_2_gene = []
    splitline2 = []
    splitline3 = []
    index = 0

    if requested_symbol in index_dict:
        requested_index = index_dict.get(requested_symbol)
    else:
        raise Exception("Incorrect requested symbol available options are:"
                        "gene_id gene_symbol gene_name entrez_id chromosome")

    # check if the value at the indicated index is over the cutoff value,
    # if it is, append the probe id
    for lines in probes_file_content:
        if lines.startswith(tuple(identifier_values1)):
            lines = lines.strip()
            splitline = lines.split(",")
            if len(splitline) > 7:
                # if the list is too long, that means there was a comma in the gene name,
                # combine them again, this won't work if there is more than 1 comma in the gene name.
                # if you get a key error, edit this to handle more commas.
                block1 = splitline[0:4]
                block2 = splitline[4:6]
                block2.insert(1, ",")
                block3 = splitline[6:]
                block2 = list(["".join(block2)])
                blocks = [block1, block2, block3]
                splitline = []

                # append them into one list again
                for lists in blocks:
                    for i in lists:
                        splitline.append(i)
            # replace the empty spots with the gene name
            # (example: when the entrez id or chromosome is missing)
            for strings in splitline:
                if not strings:
                    splitline3.append("nan" + str(index))
                else:
                    splitline3.append(strings)

            identifier_1_gene.append(splitline3[int(requested_index)].strip())
            splitline3 = []

        # repeat for the second given identifier
        if lines.startswith(tuple(identifier_values2)):
            lines = lines.strip()
            splitline = lines.split(",")
            if len(splitline) > 7:
                # if the list is too long, that means there was a comma in the gene name,
                # combine them again, this won't work if there is more than 1 comma in the gene name.
                # if you get a key error, edit this to handle more commas.
                block1 = splitline[0:4]
                block2 = splitline[4:6]
                block2.insert(1, ",")
                block3 = splitline[6:]
                block2 = list(["".join(block2)])
                blocks = [block1, block2, block3]
                splitline = []

                # append them into one list again
                for lists in blocks:
                    for i in lists:
                        splitline.append(i)
            # replace the empty spots with the gene name
            # (example: when the entrez id or chromosome is missing)
            for strings in splitline:
                if not strings:
                    splitline2.append("nan" + str(index))
                else:
                    splitline2.append(strings)
            index += 1
            identifier_2_gene.append(splitline2[int(requested_index)].strip())
            splitline2 = []

    return set(identifier_1_gene), set(identifier_2_gene)


def unique_or_shared(identifier_1_gene, identifier_2_gene):
    """
    take the extracted genes, and create new lists with the shared genes, and unique genes
    :param identifier_1_gene: list of genes for identifier 1
    :param identifier_2_gene: list of genes for identifier 2
    :return: list of shared genes, list of unique genes 1 and 2
    """

    # prepare lists for use
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

    # get the argument and set the variable
    identifier = [arguments[2], arguments[3]]
    cutoff_value = arguments[4]
    requested_symbol = arguments[1]

    # print lines to track progress
    print("getting identifiers...")
    lines_identifiers = extract_sample_annot(identifier)
    print("extracting genes...")
    gene_probes_dictionary = extract_genes()
    print("calculating highest value probes...")

    # check if there already is a pre-made list of high value probes.
    if os.path.isfile("High_Value_Probes.csv"):
        print("pre existing file found, using already improved probes")
        high_value_probes = use_premade_probes()
    else:
        print("no pre existing high value probes file found, creating one, this can take a while.")
        microarray_content = read_microarray_content()
        high_value_probes = extract_microarray_content(gene_probes_dictionary, microarray_content)
        write_intermediate_file(high_value_probes)

    print("removing values below cutoff...")
    identifier_values1, identifier_values2 = \
        extract_above_cutoff(lines_identifiers, high_value_probes, cutoff_value)

    print("extracting gene names...")
    probes_file_content = read_probes_content()
    identifier_1_gene, identifier_2_gene = \
        gene_name_finder(identifier_values1, identifier_values2,
                         requested_symbol, probes_file_content)

    print("checking similarities between genes...")
    common_genes, unique_genes_1, unique_genes_2 = \
        unique_or_shared(identifier_1_gene, identifier_2_gene)
    print("done!")

    print(
        f'{len(unique_genes_1)} unique {requested_symbol} for {arguments[2]}: {unique_genes_1} '
        f'\n{len(unique_genes_2)} unique {requested_symbol} for {arguments[3]}: {unique_genes_2}')
    print(len(common_genes), "shared", requested_symbol, common_genes)


# to protect against problems when imported
if __name__ == "__main__":
    main()
