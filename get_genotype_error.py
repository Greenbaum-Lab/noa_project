from pipeline_runner import compute_genotype_error
from utils import args_parser

if __name__ == '__main__':
    arguments = args_parser()
    compute_genotype_error(arguments.repeated, arguments.input, arguments.output)