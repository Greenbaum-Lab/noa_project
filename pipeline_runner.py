import argparse
import itertools
import numpy as np
import os
import pandas as pd
from utils import read_df_file


def mac_matrix2similarity_matrix(input_matrix):
    """
    Get numpy array of input (count of ones of the alleles in 012 format). row per individual, column per site.
    Return similarity matrix (row for individual, column for individual)
    """
    input_matrix[input_matrix == -1] = np.nan
    is_valid_window = (~np.isnan(input_matrix)).astype(np.uint8)

    window_pairwise_counts = is_valid_window @ is_valid_window.T
    # Trick to avoid computation of nan values one by one
    window0 = input_matrix.copy()
    window0[np.isnan(window0)] = 0
    window2 = input_matrix.copy()
    window2[np.isnan(window2)] = 2

    # Compute similarity
    # first_element = (ref_freq * window0) @ window0.T
    first_element = window0 @ window0.T
    # second_element = (non_ref_freq * (2 - window2)) @ (2 - window2).T
    second_element = (2 - window2) @ (2 - window2).T
    similarity = (first_element + second_element) / 4
    similarity = similarity / window_pairwise_counts
    np.fill_diagonal(similarity, -1)
    return similarity


def similarity_to_csv(similarity_matrix, individual_names, output_path):
    df = pd.DataFrame(data=similarity_matrix, columns=individual_names, index=individual_names)
    df.to_csv(output_path, index_label="Name")


def structure2genotype(file_path):
    if not os.path.exists(file_path):
        if not os.path.isdir(os.path.dirname(file_path)):
            print(f"Probably full path was not given. There is not directory (folder) name {os.path.dirname(file_path)}")
        assert False, f"File {file_path} is missing! Did you use a full path?"

    with open(file_path, "rb") as f:
        text = f.readlines()
    text = [t.decode() for t in text]
    individual_names = [t.split()[0] for t in text[1:]]
    text = [t.split() for t in text]
    assert all([len(text[0]) == len(text[i]) - 2 for i in range(1, len(text))]), "ERROR Parsing STRUCTURE file"
    text = text[1:]
    assert len(text) % 2 == 0, "ERROR Parsing STRUCTURE file"
    genotypes = [t[2:] for t in text]
    genotypes = np.array(genotypes, dtype=int)
    return genotypes, individual_names


def genotype2_mac_matrix(genotypes, individual_names):
    output_individual_names = []
    matrix012 = np.empty(shape=(genotypes.shape[0]//2, genotypes.shape[1]), dtype=np.int8)
    for col in range(genotypes.shape[1]):
        uniqs = np.unique(genotypes[:, col])
        pos_uniqs = uniqs > 0
        assert np.sum(pos_uniqs) == 2, "Not a Biallelic Data! "
    for i in range(0, genotypes.shape[0], 2):
        assert individual_names[i].startswith('"'), "ERROR Parsing STRUCTURE file"
        assert individual_names[i+1].endswith('"'), "ERROR Parsing STRUCTURE file"
        assert individual_names[i][1:] == individual_names[i+1][:-1], "ERROR Parsing STRUCTURE file"
        output_individual_names.append(individual_names[i][1:])
        count = genotypes[i: i+2] == genotypes[0]
        matrix012[i//2] = count.sum(axis=0)
        matrix012[i//2][genotypes[i: i+2].min(axis=0) < 0] = -1
    assert np.min(matrix012) >= -1, "ERROR Parsing STRUCTURE file"
    assert np.max(matrix012) <= 2, "ERROR Parsing STRUCTURE file"
    return matrix012.astype(float), output_individual_names

def filter_bad_snps_and_samples(file_path, output_dir):
    raw_data = read_df_file(file_path)
    num_of_indv = len(raw_data.ID)
    num_of_snps = len(raw_data.columns) - 1  # -1 for ID column
    snps_to_reomve = []
    for c in raw_data.columns:
        snp_fails_rate = sum(raw_data[c] == 'FAIL') / num_of_indv
        if snp_fails_rate >= 0.34:
            snps_to_reomve.append(c)
    individuals_to_remove = []
    for _, row in raw_data.iterrows():
        fails_rate = list(row).count('FAIL') / num_of_snps
        if fails_rate >= 0.34:
            individuals_to_remove.append(row.ID)
    os.makedirs(output_dir, exist_ok=True)
    with open(output_dir + 'individual_removed_by_fails.txt', "w+") as f:
        for indv in individuals_to_remove:
            f.write(f'{indv}\n')
    with open(output_dir + 'snps_removed_by_fails.txt', "w+") as f:
        for snp in snps_to_reomve:
            f.write(f'{snp}\n')
    # filtered_indv_df = raw_data[~raw_data["ID"].isin(individuals_to_remove)]
    # filtered_df = filtered_indv_df.loc[:, ~filtered_indv_df.columns.isin(snps_to_reomve)]
    filtered_df = raw_data.loc[:, ~raw_data.columns.isin(snps_to_reomve)]
    filtered_path = output_dir + 'filtered_by_bad_snps.csv'
    filtered_df.to_csv(filtered_path, index=False)
    return filtered_path

def remove_samples_by_fails(file_path, output):
    unfiltered_data = read_df_file(file_path)
    with open(output + 'individual_removed_by_fails.txt', 'r') as f:
        individuals_to_remove_list = f.read().split('\n')
    filtered_df = unfiltered_data[~unfiltered_data["ID"].isin(individuals_to_remove_list)]
    filtered_path = output + 'filtered_by_bad_samples.csv'
    filtered_df.to_csv(filtered_path, index=False)




def assign_ref_allele(data_df):
    snps_names = list(data_df.columns)
    snps_names.remove("ID")
    name_to_ref = {}
    for snp_name in snps_names:
        snp_data = data_df[snp_name]
        for i in range(len(snp_data)):
            if snp_data[i] == 'FAIL':
                continue
            else:
                name_to_ref[snp_name] = snp_data[i][0]
                break
    return name_to_ref

def genepop2012matrix(df, name_to_ref):
    minor_allele_count_df = df.copy()
    snps_names = list(df.columns)
    snps_names.remove("ID")
    for snp_name in snps_names:
        if snp_name in name_to_ref.keys():
            minor_allele_count_df[snp_name] = minor_allele_count_df[snp_name].apply(lambda x: x.count(name_to_ref[snp_name]) if x != "FAIL" else np.nan)
        else:
            minor_allele_count_df[snp_name] = minor_allele_count_df[snp_name].apply(lambda x: np.nan)
    return minor_allele_count_df

def compute_genotype_error(reapeted_file, input, output):
    """
    Compute the genotype error value. There are 2 formulas
    Args:
        reapeted_file:
        input:
        output:

    Returns:

    """
    repeat = read_df_file(reapeted_file).dropna()
    if "ID" in list(repeat.iloc[0]):
        repeat = repeat.drop([min(repeat.index)])
    data = pd.read_csv(input)
    name_to_ref = assign_ref_allele(data)
    df_per_rep = {}
    for rep in repeat.columns:
        rep_indv_names = list(repeat[rep])
        rep_data = data[data["ID"].isin(list(rep_indv_names))]
        mac_matrix = genepop2012matrix(rep_data, name_to_ref)
        df_per_rep[rep] = mac_matrix.reset_index()
    genotype_fails = {}
    snp_names = list(data.columns)
    snp_names.remove('ID')
    combinations_array = np.array(list(itertools.combinations(np.arange(len(repeat.columns)), 2)))
    for snp in snp_names:
        snp_df = pd.DataFrame(data={rep_name: df[snp] for rep_name, df in df_per_rep.items()})
        snp_matrix = snp_df.to_numpy()
        comparisons_matrix = np.abs(snp_matrix[:, combinations_array[:, 0]] - snp_matrix[:, combinations_array[:, 1]]) / 2
        genotype_fails[snp] = [np.nanmean(comparisons_matrix)]

    fails_output = output + 'genotype_errors.csv'
    pd.DataFrame(genotype_fails).to_csv(fails_output)


def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", help="input file in genepop (xlsx) format with raw sequencing")
    parser.add_argument("-o", "--output", dest="output", help="NAme of output directory. Program will generate a new "
                                                              "directory with that name, and all output files will be"
                                                              " there")
    parser.add_argument("-r", "--repeated", dest="repeated",  help="xlsx file of the repetition sequencing"
                                                                   " (only individuals name)")
    parser.add_argument("--args", dest="args", help="Any additional args")
    options = parser.parse_args()
    return options


if __name__ == '__main__':
    arguments = args_parser()
    filtered_path = filter_bad_snps_and_samples(arguments.input, arguments.output)
    compute_genotype_error(arguments.repeated, filtered_path, arguments.output)
    remove_samples_by_fails(arguments.output + 'filtered_by_bad_snps.csv', arguments.output)
    # matrix012, individual_names = genotype2_mac_matrix(genotypes, individual_names)
    # similarity_matrix = mac_matrix2similarity_matrix(matrix012)
    # similarity_to_csv(similarity_matrix, individual_names, output_file_path)
    # print("Done!")