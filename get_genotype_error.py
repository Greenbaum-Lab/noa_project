from pipeline_runner import compute_genotype_error
from utils import args_parser
import os

if __name__ == '__main__':
    arguments = args_parser()
    os.makedirs(arguments.output, exist_ok=True)
    compute_genotype_error(arguments.repeated, arguments.input, arguments.output)

    print("Done!")