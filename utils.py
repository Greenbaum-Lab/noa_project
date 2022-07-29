import argparse

import pandas as pd

def read_df_file(f_path):
    if f_path.endswith(".xlsx"):
        data = pd.read_excel(f_path)
    elif f_path.endswith(".csv"):
        data = pd.read_csv(f_path)
    else:
        assert False, "ERROR Parsing GENEpop format file"
    return data


def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", help="input file in genepop (xlsx) format with raw sequencing")
    parser.add_argument("-o", "--output", dest="output", help="Name of output directory. Program will generate a new "
                                                              "directory with that name, and all output files will be"
                                                              " there")
    parser.add_argument("-r", "--repeated", dest="repeated",  help="xlsx file of the repetition sequencing"
                                                                   " (only individuals name)")
    parser.add_argument("--args", dest="args", help="Any additional args")
    options = parser.parse_args()
    if options.output[-1] != '/':
        options.output += '/'
    return options
