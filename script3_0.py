"""
Takes information from a multitude of files, and returns the information the user requested
Author: Delshad Vegter
date: 01/11/2022
version: 5.5

To use, enter into the commandline:
            python3 script3_0.py (wanted_symbol) identifier1 identifier2 cutoff_value
example: python3 script3_0.py gene_symbol LHM PHA 17

possible requested symbols: gene_symbol, gene_id, gene_name, entrez_id, chromosome
if the found gene does not have a chromosome or entrez id,-
-they get a replacement nan value to still show comparison.
"""
import sys  # for commandline arguments
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


def extract_microarray_content(gene_probes_dictionary):
    """
    Extract the highest value probes from the microarray for use
    :param gene_probes_dictionary: dict of lists
    :return: list of lists
    """

    # predefined variables
    microarray_data = {}
    highest_value_probe_gene = []
    high_value_probes = []
    row_values = []

    # get the file content
    microarray_content = read_microarray_content()

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
        best_probe_index = row_values.index(max(row_values))
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


def use_pre_made_probes():
    """
    if a pre-made list with high value probes already exists,
    use that list instead of creating a new one again.
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
    list_identifiers = [lines_identifier1, lines_identifier2]

    with open("./data/SampleAnnot.csv", "r") as file:
        sample_annot = list(file.readlines())

    for index in range(len(identifier)):
        count = 0
        for line in sample_annot[1:]:
            count += 1
            if identifier[index] in line:
                list_identifiers[index].append(count)

    line_identifiers = [lines_identifier1, lines_identifier2]
    # check if there is anything in the list
    # if the list is empty, an invalid structure acronym was given
    for lists in line_identifiers:
        index_error = line_identifiers.index(lists)
        if not lists:
            raise Exception("no matches found for ", identifier[index_error],
                            "was there a typo? (program is case sensitive)")
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
    :param line_identifiers: list of strings
    :param high_value_probes: list of strings
    :param cutoff_value: integer
    :return: 2 lists
    """

    # predefined variables
    identifier_values1 = []
    identifier_values2 = []

    # can probably do this faster
    # go through the list and append only the values above the cutoff
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


def gene_name_finder(identifier_values1, identifier_values2, requested_symbol):
    """
    Finds the requested genes, and makes them into a list
    :param identifier_values1: list of strings
    :param identifier_values2: list of strings
    :param requested_symbol: string
    :return: 2 lists of strings
    """

    identifier_1_gene = []
    identifier_2_gene = []
    identifier_shared_genes = []
    identifier_genes = [identifier_1_gene, identifier_2_gene, identifier_shared_genes]
    index_dict = {
        "gene_id": 2,
        "gene_symbol": 3,
        "gene_name": 4,
        "entrez_id": 5,
        "chromosome": 6
    }
    splitline3 = []
    index = -1

    probes_file_content = read_probes_content()

    if requested_symbol in index_dict:
        requested_index = index_dict.get(requested_symbol)
    else:
        raise Exception("Incorrect requested symbol available options are:"
                        "gene_id gene_symbol gene_name entrez_id chromosome")

    # check if the value at the indicated index is over the cutoff value,
    # if it is, append the probe id
    for lines in probes_file_content[1:]:
        index += 1
        # check which list it should be appended to
        # so shared, id1 or id2
        # if it isn't in any of the lists, the value is set to 3,
        # and it repeats until it is in one of the lists
        if lines.startswith(tuple(identifier_values1)) and \
                lines.startswith(tuple(identifier_values2)):
            correct_list = 2
        elif lines.startswith(tuple(identifier_values1)):
            correct_list = 0
        elif lines.startswith(tuple(identifier_values2)):
            correct_list = 1
        else:
            correct_list = 3

        if correct_list < 3:
            lines = lines.strip()
            splitline = lines.split(",")
            if len(splitline) > 7:
                # if the list is too long, that means there was a comma in the gene name
                # this combines them again
                block1 = splitline[0:4]
                block2 = splitline[4:-2]
                block2.insert(1, ",")
                block3 = splitline[-2:]
                block2 = list(["".join(block2)])
                blocks = [block1, block2, block3]
                splitline = []
                # append them into one list again
                for lists in blocks:
                    for i in lists:
                        splitline.append(i)
            # replace the empty spots with nan+index
            # (example: when the entrez id or chromosome is missing)
            for strings in splitline:
                if not strings:
                    splitline3.append("nan" + str(index))
                else:
                    splitline3.append(strings)
            identifier_genes[correct_list].append(splitline3[int(requested_index)].strip())
            splitline3 = []

    return set(identifier_1_gene), set(identifier_2_gene), set(identifier_shared_genes)


def main():
    """
    Main, where all functions are run in order
    :return: /
    """

    arguments = sys.argv

    if len(arguments) < 5:
        raise Exception(
            "too few arguments, please use: "
            "(python3 script2_3.py wanted_symbol identifier1 identifier2 cutoff_value)")
    try:
        # if this runs it's a valid input
        cutoff_value = float(arguments[4])
    except ValueError as invalid_cutoff_value:
        # if it gives a ValueError, report that to the user
        raise Exception("invalid cutoff value") from invalid_cutoff_value

    # get the argument and set the variable
    identifier = [arguments[2], arguments[3]]
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
        high_value_probes = use_pre_made_probes()
    else:
        print("no pre existing high value probes file found, creating one, this can take a while.")
        high_value_probes = extract_microarray_content(gene_probes_dictionary)
        write_intermediate_file(high_value_probes)

    print("removing values below cutoff...")
    identifier_values1, identifier_values2 = \
        extract_above_cutoff(lines_identifiers, high_value_probes, cutoff_value)

    print("extracting gene names...")
    identifier_1_gene, identifier_2_gene, identifier_shared_genes = \
        gene_name_finder(identifier_values1, identifier_values2, requested_symbol)
    print("done!")

    print(
        f'{len(list(identifier_1_gene))} unique {requested_symbol} '
        f'for {arguments[2]}: {list(identifier_1_gene)} '
        f'\n{len(list(identifier_2_gene))} unique {requested_symbol} '
        f'for {arguments[3]}: {identifier_2_gene}')
    print(len(identifier_shared_genes), "shared", requested_symbol, identifier_shared_genes)


# to protect against problems when imported
if __name__ == "__main__":
    main()
