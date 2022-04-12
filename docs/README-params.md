
# Parameters

## Moltraits

- `MOLTRAITS_ID_COL`: Column in molecular trait input file (`moltraits.tsv.gz`) containing feature IDs.
- `MOLTRAITS_SAMPLE_START`: Column in molecular trait input file (`moltraits.tsv.gz`) corresponding to index where sample data starts. Index should be 0-indexed. All feature metadata should be before sample start.

## Covariates

- `COVARIATES`: List of columns in covariates to include in models.
- `COVARIATES_REFERENCES`: List of reference values for categorical covariates. Syntax: ["cov_1::ref", "cov_2::ref", "cov_3::ref"]
- `COVARIATES_MEAN_IMPUTE`: List of columns in covariates to mean impute if missing.

## Covariance matrices

- `MATRICES`: Dictionary of effects to test and corresponding data files. The key should be the name, and the value should be the path to the file. Files can either be (1) TSVs where the first column is the sample ID and the other columns are values to compute the covariance matrix with, or (2) TSVs where it is a pre-computed covariance matrix (without row names or columns). If (2), the order of the rows/columns must **match** the order of samples in the `samples.txt` input file.
- `MATRICES_TRANSFORM_DATAFRAMES`: List of effects to compute covariance matrix.

## PEER

- `PEER_N_FACTORS`: List of the number of PEER factors to adjust the expression data for. If "NA", the pipeline will not run PEER and use the data in `moltraits.tsv.gz` as direct input into LIMIX. If 0, the pipeline will still run PEER, which means if `PEER_USE_RESIDUALS` is True, we will correct for covariates using PEER instead of including covariates directly in LIMIX.
- `PEER_ITERATIONS`: Number of iterations to run PEER. Recommended: 1000
- `PEER_ACCOUNT_MEAN`: ["True" | "False"]. Add additional covariate to account for mean moltrait signal. **Note**: Value is a string.
- `PEER_INVERSE_NORMALIZE`: ["True" | "False"]. Inverse normalize the moltrait data prior to running PEER. **Note**: Value is a string.
- `PEER_USE_RESIDUALS`: [true | false]. Use PEER residuals after covariate and factor correction as input into LIMIX.

## Variance decomposition

- `VARDEC_NORMALIZATION`: Normalize moltrait data before modeling with LIMIX. Options: "NA" (no normalization), "rank_normal", 'sklearn_box_cox', 'sklearn_rank_normal' (inverse rank normalization as implemented in `sklearn`)
- `VARDEC_RESIDUAL_DISTRIBUTION`: Expected distribution of residual effect sizes in the fitted linear mixed model. Options: "normal", "bernoulli", "probit", "binomial", "poisson".
- `VARDEC_MATRICES`: Dictionary of effects to test. Each key should correspond to a label for test and values should be a list of effects to test.
