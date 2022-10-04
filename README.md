# noa_project

###Overall
Written and supported by Shahar Mazie from HUJI University (shahar.mazie@mail.huji.ac.il)

This script is filtering DNA sequenced data, and computing some pairwise metrics on it.

The inputs for the script are:

1. -o Output directory name. It can exits or not exists. The script WILL overwrite files! 
This will make your life easier if you have mid-results, since it can skip steps in case it failed in the middle.

2. -i input file name of the sequenced data. Should be csv or xlsx file. The format was provided by Noa Kan.

3. -r Repeats file [Optional]. Csv or Xlsx file by the format that was given by Noa Kan.

The script should be able to run easily with any OS and any Python 3 with some standard packages. It was tested with Python 3.8.



###Steps

1. Filter SNPs with >34% of missing data
2. Compute number of mistakes and genotype error for every SNP (using Repeats file)
3. Estimate lambda1 and lambda2 of a mixture poisson distribution. 
4. Compute threshold for errors.
5. Filter SNPs by the threshold, and compute a new genotype error based on the SNPS left with a single poisson distribution.
6. Filter samples with >34% missing data.
7. Assign a "ref" allele for every SNP and generate 012 matrix
8. Compute f_mu_matrix - EQ 4 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1461426/ 
9. Compute weighted f_mu score - Gili's score from https://academic.oup.com/genetics/article/202/4/1299/5930174?login=false

* If Repeats file is missing, we skip steps 2-5 and use the output from step 1 from step 6. 
### Output files 
All outputs will be find in the output directory provided with -o.

The output files of the script - relative to the step numbers:
1. snps_removed_by_fails.txt - names of the SNPS that were removed.
filtered_by_bad_snps.csv - csv of the data after removing the filtered SNPs

2.genotype_num_of_errors.csv - per SNP the number of errors it had.
genotype_errors.csv - per SNP its estimated genotype error.

3. mixture_poisson_params.txt - parameters of the mixture poisson (lambda1, lambda2, m)

4. snps_removed_by_genotype_error.txt - list of SNP names that were removed thanks to high genotype error

5. filtered_all_snps.csv - Data frame after filtering all SNPs

6. filtered.csv - Data frame after all filters (including samples with >34% missing data).
individual_removed_by_fails.txt - names of individuals removed because of >34% missing data.

7.f_mu_matrix.csv - matrix of distances based on f_mu metric.
f_mu_matrix_pairs.csv - same matrix but each line is pair of individual and their similarity 

8.weighted_f_mu_matrix.csv - matrix of distances based on weighted f_mu metric.
weighted_f_mu_matrix_pairs.csv - same matrix but each line is pair of individual and their similarity 

###Computing Similarity

It is possible to filter the data yourself and use this script just to compute a similarity matrix of the data.

In order to do that, run "compute_similarity.py" script (instead of pipeline_runner.py, because, well, 
we don't run the whole pipeline).
Give with -i the path to your data csv file like in the pipeline_runner, and -o with output path of an directory (might be not exists yet).
The script will write the files in points 7,8 in "Output files" section in your output directory. 
