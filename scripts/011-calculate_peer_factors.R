#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c("-i", "--phenotypes"),
                        type="character",
                        default="phenotypes.tsv.gz",
                        help="
                            TSV / BED file containing phenotype information.
                            Must have the following
                            column names:
                            chr
                            start
                            end
                            <phenotype-id>: specified with `p-id` flag
                            <sample_names>: columns for sample names
                        "
  ),

  optparse::make_option(c("-l", "--sample_column_start"),
                        type="numeric",
                        default=7,
                        help="
                            Column in which the sample ids begin (1 based).
                            For instance(chr, start, end, feature_id, sample_1,
                            sample_2, ...) would be 5, and (feature_id,
                            sample_1, sample_2, ...) would be 2.
                            [default %default]
                        "
  ),

  optparse::make_option(c("-k", "--hidden_factors"),
                        type="numeric",
                        default=15,
                        help="The number of hidden factors. [default %default]"
  ),

  optparse::make_option(c("-n", "--iterations"),
                        type="numeric",
                        default=1000,
                        help="Maximum number of iterations. [default %default]"
  ),

  optparse::make_option(c("-c", "--covariate_file"),
                        type="character",
                        default=NULL,
                        help="
                            Input covariates file. This should be a matrix
                            where samples are rows and covariates are columns.
                            [default %default]
                        "
  ),

  optparse::make_option(c("-m", "--account_mean"),
                        type="logical",
                        default=TRUE,
                        help="Add additional covariate to account for mean y
                        signal. [default %default]"
  ),

  optparse::make_option(c("-z", "--inverse_norm"),
                        type="logical",
                        default=TRUE,
                        help="
                            Inverse normalize y features before running PEER
                            (so each y feature has a normal distribution across
                            samples). Note if true, you do not need to include
                            a mean covariate. [default %default]
                        "
  ),

  optparse::make_option(c("-o", "--output_base"),
                        type="character",
                        default="moltraits",
                        help="Prefix for output files. [default %default]"
  ),

  optparse::make_option(c("-s", "--seed"),
                        type="numeric",
                        default=1938,
                        help="Random number seed. [default %default]"
  ),

  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = TRUE,
                        help = ""
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Corrects diffential results using IHW."
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
  optionStrings <- character()
  for (item in parserObj@options) {
    optionStrings <- append(optionStrings,
                            c(item@short_flag, item@long_flag))
  }
  optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(peer))
suppressPackageStartupMessages(library(dplyr))
################################################################################

################################ Functions #####################################
invnorm_random = function(x) {
  qnorm((rank(x, na.last="keep", ties.method="random") - 0.5) / sum(!is.na(x)))
}

calculate_peer_factors <- function(mtrx,
                                   k = 10,
                                   iterations,
                                   output_base,
                                   other_dat=NULL,
                                   covariates_mtx=NULL,
                                   include_mean=T,
                                   inverse_norm_before=T,
                                   verbose = T) {
  # Input:
  # mtrx - matrix of samples (rows) by phenotypes (cols)
  # other_dat - dataframe of samples (cols) by feature_ids (rows)
  # covariates_mtx - matrix of samples (rows) by feature_ids (cols)

  if (inverse_norm_before) {
    cat("Inverse normalizing each feature id (col) across samples (row).\n")
    # apply invnorm_random to each gene (column) across samples (rows)
    mtrx = as.matrix(apply(mtrx, 2, invnorm_random))
    cat("\nstandard deviation cols:\n")
    print(head(apply(mtrx, 2, sd)))
    cat("\nstandard deviation rows:\n")
    print(head(apply(mtrx, 1, sd)))
    cat("\n")
  }
  print(paste("PEER input matrix shape:", paste(dim(mtrx), collapse=",") ))

  # init peer
  model = peer::PEER()

  # Set observed data
  peer::PEER_setPhenoMean(model, mtrx)

  # set number of hidden factors
  peer::PEER_setNk(model, k)
  if (k > 0) {
    factor_names = paste("factor", seq(1,k), sep="_")
    factor_names_full = factor_names
  } else {
    factor_names = c()
    factor_names_full = c()
  }
  print(paste("PEER_factors=", (k), sep=""))

  # set number of max iterations
  PEER_setNmax_iterations(model, iterations)

  # For most use cases, including the mean effect is likely to be a good choice.
  # If y is inverse normalized including a mean effect does not make sense.
  if (include_mean) {
    factor_names_full = c("mean", factor_names_full)
  }
  PEER_setAdd_mean(model, include_mean)

  # Measured experimental variables that may contribute to variability in the
  # data can be included in the inference as covariates
  # (e.g batch, RNA quality, ...)
  if (!is.null(covariates_mtx)) {
    cat("Setting covariates:", colnames(covariates_mtx), "\n")
    factor_names_full = c(colnames(covariates_mtx), factor_names_full)
    #print(factor_names_full)
    # if there is an error it is probly because covariates_mtx is not a matrix
    PEER_setCovariates(model, covariates_mtx)
  }

  # Perform inference > results in a model object with posterior distributions
  # of the variables
  # aka now run peer
  PEER_update(model)

  # get the posterior mean of the inferred factors (NxK matrix),
  # their weights (GxK matrix),
  # precision (inverse variance) of the weights (Kx1 matrix),
  # the residual dataset (NxG matrix):
  factors = PEER_getX(model)
  weights = PEER_getW(model)
  precision = PEER_getAlpha(model)
  residuals = PEER_getResiduals(model)

  # add labels
  #print(factor_names_full) # print
  feature_ids = colnames(mtrx)
  #samples = unlist(lapply(rownames(mtrx), FUN=function(x) {
  #    strsplit(x, "X")[[1]][2] } ))
  samples = rownames(mtrx)

  factors = data.frame(t(factors))
  colnames(factors) = samples
  rownames(factors) = factor_names_full # this may give a bug
  factors$factor = factor_names_full
  # move factor column first
  factors = factors[,c(ncol(factors),1:(ncol(factors)-1))]
  cat("processed factors\n")

  precision_old = precision
  precision = data.frame(precision)
  rownames(precision) = factor_names_full
  colnames(precision) = "precision"
  precision$factor = factor_names_full
  # move factor column first
  precision = precision[,c(ncol(precision),1:(ncol(precision)-1))]
  cat("processed precision\n")

  weights = data.frame(weights)
  rownames(weights) = feature_ids
  colnames(weights) = factor_names_full
  if (!is.null(other_dat)) {
    weights = merge(other_dat, weights, by="row.names")
    weights[1] = NULL
  } else {
    weights$feature_id = feature_ids
    # move feature_ids column first
    weights = weights[,c(ncol(weights),1:(ncol(weights)-1))]
  }
  cat("processed weights\n")

  residuals = data.frame(t(residuals), check.names=FALSE)
  colnames(residuals) = samples
  rownames(residuals) = feature_ids
  # inverse normalize residuals
  residuals_inv = data.frame(t(apply(as.matrix(residuals), 1, invnorm_random)),
                             check.names=FALSE)
  # if other dat, it is assumed the colnames matched the original input mol
  # trait matrix ... (colnames given to function)
  if (!is.null(other_dat)) {
    residuals = merge(other_dat, residuals, by="row.names")
    residuals[1] = NULL
    residuals_inv = merge(other_dat, residuals_inv, by="row.names")
    residuals_inv[1] = NULL
  } else {
    # else, add a column called feature id and give it the rownames of the
    # input matrix
    residuals$feature_id = feature_ids
    # move probe column first
    residuals = residuals[,c(ncol(residuals),1:(ncol(residuals)-1))]

    residuals_inv$feature_id = feature_ids
    residuals_inv = residuals_inv[,c(ncol(residuals_inv),
                                     1:(ncol(residuals_inv)-1))]
    colnames(residuals_inv) = colnames(residuals)
  }
  cat("processed residuals\n")

  # save our matricies
  sep="\t"
  gzfh = gzfile(paste(output_base, "-peer_factors.tsv.gz", sep = ""), "w")
  write.table(factors, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)

  gzfh = gzfile(paste(output_base, "-peer_weights.tsv.gz", sep = ""), "w")
  write.table(weights, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)

  gzfh = gzfile(paste(output_base, "-peer_precision.tsv.gz", sep = ""), "w")
  write.table(precision, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)

  gzfh = gzfile(paste(output_base, "-peer_residuals.tsv.gz", sep = ""), "w")
  write.table(residuals, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)

  gzfh = gzfile(paste(output_base, "-peer_residuals-invnorm.tsv.gz", sep = ""),
                "w")
  write.table(residuals_inv, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
set.seed(arguments$options$seed)

# Get arguments
if (verbose) {
  print("Reading in the data...")
}

phenotypes <- read.csv(
  arguments$options$phenotypes,
  sep = "\t",
  header = T
)
# subset columns
smpl_start <- arguments$options$sample_column_start
# rownames(dat) <- paste("row", seq(1,nrow(dat)), sep="_")
pheno_mtx <- phenotypes[,colnames(phenotypes)[smpl_start:ncol(phenotypes)]]
pheno_mtx <- t(pheno_mtx)
pheno_other <- phenotypes[,colnames(phenotypes)[1:smpl_start-1]]

# Get covariates
if (!is.null(arguments$options$covariate_file)) {
  covariates_mtx <- read.csv(
    arguments$options$covariate_file,
    sep = "\t",
    header = T,
    row.names = 1
  )
  if (!all(rownames(covariates_mtx) == rownames(pheno_mtx))) {
    stop(paste(
        "The samples (order) in the covariates matrix don't match the",
        "samples (order) in the phenotype matrix."
    ))
  }
  covariates_mtx <- as.matrix(covariates_mtx)
} else {
  covariates_mtx <- NULL
}

# run peer
calculate_peer_factors(
  pheno_mtx,
  arguments$options$hidden_factors,
  arguments$options$iterations,
  arguments$options$output_base,
  pheno_other,
  covariates_mtx,
  arguments$options$account_mean,
  arguments$options$inverse_norm
)

if (verbose) {
  print("Done.")
}

################################################################################
