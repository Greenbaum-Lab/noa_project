import warnings

import numpy as np
import os
import pandas as pd
from tqdm import tqdm

from utils import read_df_file, args_parser
import matplotlib.pyplot as plt
from scipy.special import factorial


def mac_matrix2similarity_matrix(input_matrix, weighted=False):
    """
    Get numpy array of input (count of ones of the alleles in 012 format). row per individual, column per site.
    Return similarity matrix (row for individual, column for individual)
    """
    input_matrix[input_matrix == -1] = np.nan
    is_valid_window = (~np.isnan(input_matrix)).astype(np.uint16)
    window_pairwise_counts = is_valid_window @ is_valid_window.T
    # Trick to avoid computation of nan values one by one
    window0 = input_matrix.copy()
    window0[np.isnan(window0)] = 0
    window2 = input_matrix.copy()
    window2[np.isnan(window2)] = 2

    if weighted:
        num_valid_genotypes = np.sum(is_valid_window, axis=0)
        non_ref_count = np.sum(input_matrix == 1, axis=0) + 2 * np.sum(input_matrix == 2, axis=0)
        non_ref_freq = non_ref_count / (2 * num_valid_genotypes)
        ref_freq = 1 - non_ref_freq
        first_element = (ref_freq * window0) @ window0.T
        second_element = (non_ref_freq * (2 - window2)) @ (2 - window2).T
    else:
        first_element = window0 @ window0.T * 2
        second_element = (2 - window2) @ (2 - window2).T

    similarity = (first_element + second_element) / 4
    similarity = similarity / window_pairwise_counts
    np.fill_diagonal(similarity, -1)
    return similarity


def similarity_to_csv(similarity_matrix, individual_names, output_path):
    if os.path.exists(output_path):
        return
    result_file = ""
    for idx1, name1 in enumerate(individual_names):
        for idx2, name2 in enumerate(individual_names[idx1 + 1:], start=idx1+1):
            result_file += f"{name1},{name2},{similarity_matrix[idx1, idx2]}\n"
    with open(output_path + '_pairs.csv', 'w') as f:
        f.write(result_file[:-1])
    df = pd.DataFrame(data=similarity_matrix, columns=individual_names, index=individual_names)
    df.to_csv(output_path + '.csv', index_label="Name")


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


def filter_bad_samples(output):
    filtered_path = output + 'filtered.csv'
    if os.path.exists(filtered_path):
        return filtered_path
    snps_filtered_data = read_df_file(output + 'filtered_all_snps.csv')
    num_of_snps = len(snps_filtered_data.columns) - 1  # -1 for ID column
    individuals_to_remove = {}
    for _, row in snps_filtered_data.iterrows():
        fails_rate = list(row).count('FAIL') / num_of_snps
        if fails_rate >= 0.34:
            individuals_to_remove[row.ID] = fails_rate
    with open(output + 'individual_removed_by_fails.txt', "w+") as f:
        for indv, val in individuals_to_remove.items():
            f.write(f'{indv}, {val}\n')
    filtered_df = snps_filtered_data[~snps_filtered_data["ID"].isin(set(list(individuals_to_remove.keys())))]
    filtered_df.to_csv(filtered_path, index=False)
    return filtered_path


def filter_bad_SNPs(file_path, output_dir):
    filtered_path = output_dir + 'filtered_by_bad_snps.csv'
    if os.path.exists(filtered_path):
        return filtered_path
    raw_data = read_df_file(file_path)
    num_of_indv = len(raw_data.ID)
    snps_to_reomve = {}
    for c in raw_data.columns:
        snp_fails_rate = sum(raw_data[c] == 'FAIL') / num_of_indv
        if snp_fails_rate >= 0.34:
            snps_to_reomve[c] = snp_fails_rate

    with open(output_dir + 'snps_removed_by_fails.txt', "w+") as f:
        for snp, val in snps_to_reomve.items():
            f.write(f'{snp}, {val}\n')
    filtered_df = raw_data.loc[:, ~raw_data.columns.isin(set(list(snps_to_reomve.keys())))]
    filtered_df.to_csv(filtered_path, index=False)
    return filtered_path


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

def compute_genotype_error(reapeted_file, input_name, output):
    """
    Compute the genotype error value by Alan computation.
    Args:
        reapeted_file:
        input:
        output:

    Returns:

    """
    fails_num_output = output + 'genotype_num_of_errors.csv'
    fails_output = output + 'genotype_errors.csv'
    mle_text_path = output + "mixture_poisson_params.txt"
    repeat = read_df_file(reapeted_file).dropna()
    if "ID" in list(repeat.iloc[0]):
        repeat = repeat.drop([min(repeat.index)])
    data = pd.read_csv(input_name)
    name_to_ref = assign_ref_allele(data)
    df_per_rep = {}
    num_of_repeats = len(list(repeat.columns))
    for rep in repeat.columns:
        rep_indv_names = list(repeat[rep])
        rep_data = data[data["ID"].isin(list(rep_indv_names))]
        mac_matrix = genepop2012matrix(rep_data, name_to_ref)
        df_per_rep[rep] = mac_matrix.reset_index()
    d3_mat = np.concatenate([np.expand_dims(pd.DataFrame.to_numpy(x), axis=0) for x in df_per_rep.values()], axis=0)
    d3_mat = d3_mat[:, :, 1:]
    num_of_repeated_individuals = d3_mat.shape[1]
    if all(os.path.exists(p) for p in [fails_output, fails_num_output, mle_text_path]):
        return num_of_repeated_individuals, num_of_repeats
    for i in range(d3_mat.shape[1]):
        assert np.all(d3_mat[:, i, 0] == list(repeat.iloc[i])), f"line {i} is not aligned between repeats and matrix." \
                                                                f"Developer bug. Contact Shahar"
    d3_mat = d3_mat[:, :, 1:]
    snp_names = list(data.columns)
    genotype_fails = {}
    genotype_num_fails = {}
    snp_names.remove('ID')
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in reduce")
        max_min_diff = (np.max(d3_mat, axis=0) - np.min(d3_mat, axis=0)).astype(float)
    same_counter = np.count_nonzero(max_min_diff == 0, axis=0)
    non_nan_counter = np.count_nonzero(~np.isnan(max_min_diff), axis=0)
    num_of_misses = non_nan_counter - same_counter
    mle_brute_force(num_of_misses, mle_text_path)
    for idx, snp_name in enumerate(snp_names):
        genotype_fails[snp_name] = [1 - (same_counter[idx] / non_nan_counter[idx]) ** (1 / num_of_repeats)]
        genotype_num_fails[snp_name] = [num_of_misses[idx]]

    df = pd.DataFrame.from_dict(genotype_fails)
    df.to_csv(fails_output, index=False)

    df = pd.DataFrame.from_dict(genotype_num_fails)
    df.to_csv(fails_num_output, index=False)


def compute_threshold_and_filter_by(output, num_of_repeated_individuals, num_of_repeats):
    filtered_all_snps_path = output + 'filtered_all_snps.csv'
    if os.path.exists(filtered_all_snps_path):
        return
    input_name = output + 'genotype_num_of_errors.csv'
    data = read_df_file(input_name)
    missing_in_data = np.array(data.iloc[0]).astype(int)
    MAX_MISSING = 21
    k = np.arange(MAX_MISSING + 1)
    with open(output + "mixture_poisson_params.txt", "r") as f:
        mle_txt = f.read().split()

    e1_prefer = admixture_const_prob(lambda1=float(mle_txt[0]), lambda2=float(mle_txt[1]), m=float(mle_txt[2]), k=k)
    max_num_of_accepted_errors = np.max(np.where(e1_prefer > 0)[0])
    print(f"{max_num_of_accepted_errors} is the highest value of misses to be accepted.")
    datat = data.T
    filter_SNP_names = list((datat[datat[0] > max_num_of_accepted_errors]).index)
    with open(output + "snps_removed_by_genotype_error.txt", "w") as f:
        f.write('\n'.join(filter_SNP_names))
    half_filtered_data = read_df_file(output + 'filtered_by_bad_snps.csv')
    filter_ge_snps = half_filtered_data.drop(filter_SNP_names, axis=1)
    filter_ge_snps.to_csv(filtered_all_snps_path, index=False)
    mle_for_single_poisson = np.mean(missing_in_data[missing_in_data <= max_num_of_accepted_errors])
    est_ge = 1 - ((num_of_repeated_individuals - mle_for_single_poisson) / num_of_repeated_individuals) ** (1 / num_of_repeats)
    print(f"Estimated genotype error for remaining SNPS is {est_ge}")


def mle_brute_force(data, output_path):
    current_mle = -np.inf
    mle_params = [None, None, None]
    factorial_data = factorial(data)
    lambda_1_center = 2.01
    lambda_2_center = 10
    lambda_1_eps = [2, 0.4, 0.05]
    lambda_2_eps = [9, 2, 0.1]
    print("Computing MLE for multi poisson distribution - this will take a few minutes")
    for iter in tqdm(range(3)):
        for lambda1 in tqdm(np.linspace(start=lambda_1_center - lambda_1_eps[iter], stop=lambda_1_center + lambda_1_eps[iter], num=100), leave=False):
            for lambda2 in np.linspace(start=lambda_2_center - lambda_2_eps[iter], stop=lambda_2_center + lambda_2_eps[iter], num=100):
                if np.min([lambda1, lambda2]) <= 0:
                    continue
                for m in np.linspace(start=0, stop=1, num=100 * (iter + 1)):
                    res = apply_admixture_prob_func(lambda1, lambda2, m, data, factorial_data)
                    if res > current_mle:
                        current_mle = res
                        mle_params = [lambda1, lambda2, m]
        print(mle_params)
        lambda_1_center = mle_params[0]
        lambda_2_center = mle_params[1]

    print(f"Done! with lambda1={mle_params[0]}, lambda2={mle_params[1]}, m={mle_params[2]}")
    with open(output_path, "w") as f:
        f.write('\n'.join([str(e) for e in mle_params]))

def apply_admixture_prob_func(lambda1, lambda2, m, data, fac_data):
    e1 = m * (((lambda1 ** data) * (np.e ** - lambda1)) / fac_data)
    e2 = (1-m) * (((lambda2 ** data) * (np.e ** - lambda2)) / fac_data)
    return np.sum(np.log(e1 + e2))


def admixture_const_prob(lambda1,  lambda2, m, k=None):
    e1 = const_prob(lambda1, k, m)
    e2 = const_prob(lambda2, k, m)
    return e1 > e2

def const_prob(lambda_x, k, m):
    """
    Compute P(I,K). Followed by base rules, if P(i,k) > P(j,k) then P(i|k) > P(j|k)
    :param lambda_x:
    :param k:
    :param m:
    :return:
    """
    return m * (((lambda_x ** k) * (np.e ** - lambda_x)) / factorial(k))


def filtered_data2similarity(data, output):
    name_to_ref = assign_ref_allele(data)
    matrix012 = genepop2012matrix(data, name_to_ref)
    matrix012 = pd.DataFrame.to_numpy(matrix012.drop(["ID"], axis=1))
    f_mu_matrix_path = output + 'f_mu_matrix'
    if not os.path.exists(f_mu_matrix_path):
        f_mu_similarity_matrix = mac_matrix2similarity_matrix(matrix012, weighted=False)
        similarity_to_csv(f_mu_similarity_matrix, list(data.ID), f_mu_matrix_path)
    weighted_f_mu_matrix_path = output + 'weighted_f_mu_matrix'
    if not os.path.exists(weighted_f_mu_matrix_path):
        weighted_similarity_matrix = mac_matrix2similarity_matrix(matrix012, weighted=True)
        similarity_to_csv(weighted_similarity_matrix, list(data.ID), weighted_f_mu_matrix_path)

if __name__ == '__main__':
    arguments = args_parser()
    os.makedirs(arguments.output, exist_ok=True)
    filtered_snp_path = filter_bad_SNPs(arguments.input, arguments.output)
    if arguments.repeated:
        num_of_repeated_individuals, num_of_repeats = compute_genotype_error(arguments.repeated, filtered_snp_path, arguments.output)
        compute_threshold_and_filter_by(arguments.output, num_of_repeated_individuals, num_of_repeats)
        filtered_path = filter_bad_samples(arguments.output)
    else:
        filtered_path = filtered_snp_path
    data = read_df_file(filtered_path)
    filtered_data2similarity(data, arguments.output)
    print("Done!")
