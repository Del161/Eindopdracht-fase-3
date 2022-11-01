import sys

"""
function_placeholder
Author: Delshad Vegter 
date: 01/11/2022
version: 1
"""
def extract_microarray_content(input_name):
    """
    Opening of the given file.
    :param input_name: the input name given in commandline
    :return: list with content of the file
    """

    # predefined variables
    probe_ID = ""

    microarray_data = {}

    with open(input_name, "r") as file:
        # set file lines to variable
        file_content = list(file.readlines())

    # take each line, turn the id into a key and the content into the key content.
    for lines in file_content:
        lines = lines.split(",")
        probe_ID = lines.pop(0)
        probe_content = "".join(lines)
        microarray_data[probe_ID] = probe_content

    return microarray_data



def main():
    """
    Main, where all functions are ran in order
    :return: /
    """

    arguments = sys.argv

    if len(arguments) < 2:
        raise Exception(
            "No input filename given, please use commandline arguments. "
            "(python3 script2_3.py file inputname, output filename )")
    else:
        # get the argument and set the variable
        input_name = arguments[1]

    array_file_content = extract_microarray_content(input_name)


# to protect against problems when imported
if __name__ == "__main__":
    main()