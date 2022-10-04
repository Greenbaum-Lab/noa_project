import os

from pipeline_runner import filtered_data2similarity
from utils import args_parser, read_df_file

if __name__ == '__main__':
    arguments = args_parser()
    os.makedirs(arguments.output, exist_ok=True)
    data = read_df_file(arguments.input)
    filtered_data2similarity(data, arguments.output)