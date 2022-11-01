import sys

"""
function_placeholder
Author: Delshad Vegter 
date: 01/11/2022
version: 1.1
"""
def extract_microarray_content():
    """
    Opening of the given file, extract microarray data to dictionary.
    :return: list with content of the file
    """

    # predefined variables
    probe_ID = ""
    probe_ID_list = []
    microarray_data = {}
    index = 0

    with open("./data/MicroarrayExpression.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    # take each line, turn the id into a key and the content into the key content.
    for lines in file_content:
        lines = lines.strip().split(",")
        probe_ID = lines.pop(0)
        probe_content = "".join(lines)
        probe_ID_list.append(probe_ID)
        microarray_data[probe_ID] = probe_content

        # find which line contains the identifier
        if identifier in lines:
            locations.append(index)
        index += 1

    return microarray_data, probe_ID_list

def extract_sample_annot():


    with open("./data/SampleAnnot.csv", "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    probe_ID = ""
    probes_data = {}

    # take each line, turn the id into a key and the content into the key content.
    for lines in file_content:
        lines = lines.strip().split(",")
        probe_ID = lines.pop(0)
        probe_content = "".join(lines)
        probes_data[probe_ID] = probe_content

    return probes_data


def main():
    """
    Main, where all functions are ran in order
    :return: /
    """

    arguments = sys.argv

# if len(arguments) < 2:
#     raise Exception(
#         "No input filename given, please use commandline arguments. "
#         "(python3 script2_3.py file inputname, output filename )")
# else:
#     # get the argument and set the variable
#     input_name = arguments[1]

    array_file_content, probe_ID_list = extract_microarray_content()
    probes_data = extract_sample_annot()


# to protect against problems when imported
if __name__ == "__main__":
    main()