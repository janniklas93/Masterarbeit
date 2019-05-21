library(stringr)
library(ggplot2)
library(cowplot)
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/marker_genes.R")
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/Simulation.R")

### Helper functions
evaluate = function(row){
    ind = which.max(row)
    name = names(row)[ind]
}
evaluate2 = function(row){
    sum(row) > 0
}
evaluate3 = function(row){
    amt = sum(row > 5)
    if(amt > 1)
        return(FALSE)
    else
        return(TRUE)
}

### Read in bulk data used as basis of simulation
# Alpha-cells, beta-cells, delta-cells MOUSE
alpha_beta_delta = read.table(
    "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Purified/alpha_beta_delta_mouse/GSE76017_data_geo.tsv",
    header = TRUE,
    sep = "\t"
)
alpha_bulk = alpha_beta_delta[,grep("Alpha", colnames(alpha_beta_delta))]
beta_bulk = alpha_beta_delta[,grep("Beta", colnames(alpha_beta_delta))]
delta_bulk = alpha_beta_delta[,grep("Delta", colnames(alpha_beta_delta))]
# Alpha-cells, beta-cells HUMAN
alpha_beta = read.table(
    "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Purified/alpha_beta_purified/purified_alpha_beta_bulk_tpm.tsv",
    header = TRUE,
    sep = "\t"
)
alpha_bulk = alpha_beta[,grep("Alpha", colnames(alpha_beta))]
beta_bulk = alpha_beta[,grep("Beta", colnames(alpha_beta))]

### Read in single-cell seq data used for training
# Read in meta info
meta = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Meta_information.tsv"
meta = read.table(meta, header = TRUE, sep = "\t")
rownames(meta) = meta$Name

# Read in Lawlor
lawlor_counts = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Human_differentiated_pancreatic_islet_cells_scRNA/Lawlor.tsv"
lawlor_counts = read.table(lawlor_counts, header = TRUE, sep = "\t")
colnames(lawlor_counts) = str_replace(colnames(lawlor_counts), pattern = "\\.", "_")
# Process Lawlor
lawlor_meta = meta[colnames(lawlor_counts),]
lawlor_counts = lawlor_counts[, lawlor_meta$Subtype %in% c("Alpha", "Beta", "Gamma", "Delta")]
rownames(lawlor_counts) = str_to_upper(rownames(lawlor_counts))
lawlor_types = lawlor_meta$Subtype[rownames(lawlor_meta) %in% colnames(lawlor_counts)]
lawlor_types = as.character(lawlor_types)

# Read in Segerstolpe
segerstolpe_counts = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv"
segerstolpe_counts = read.table(segerstolpe_counts, header = TRUE, sep = "\t")
colnames(segerstolpe_counts) = str_replace(colnames(segerstolpe_counts), pattern = "\\.", "_")
# Process Segerstolpe
seg_meta = meta[colnames(segerstolpe_counts),]
segerstolpe_counts = segerstolpe_counts[, seg_meta$Subtype %in% c("Alpha", "Beta", "Gamma", "Delta")]
rownames(segerstolpe_counts) = str_to_upper(rownames(segerstolpe_counts))
seg_types = seg_meta$Subtype[rownames(seg_meta) %in% colnames(segerstolpe_counts)]
seg_types = as.character(seg_types)

# Read in Baron
baron_counts = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Human_differentiated_pancreatic_islet_cells_scRNA/Baron_human.tsv"
baron_counts = read.table(baron_counts, header = TRUE, sep = "\t")
colnames(baron_counts) = str_replace(colnames(baron_counts), pattern = "\\.", "_")
# Process Baron
baron_meta = meta[colnames(baron_counts),]
baron_counts = baron_counts[,baron_meta$Subtype %in% c("Alpha", "Beta", "Gamma", "Delta")]
rownames(baron_counts) = str_to_upper(rownames(baron_counts))
baron_types = baron_meta$Subtype[rownames(baron_meta) %in% colnames(baron_counts)]
baron_types = as.character(baron_types)

### Choose current training data, process and prepare
training_data = baron_counts
training_types = baron_types
colnames(training_data) = str_replace(colnames(training_data), pattern = "\\.", "_")
colnames(training_data) = str_replace_all(colnames(training_data), pattern = "^X", "")
# Data cleaning
row_var = apply(training_data, FUN = var, MARGIN = 1)
training_data = training_data[row_var != 0,]
training_data = training_data[rowSums(training_data) >= 1,]
# Determine marker genes by differential expression
marker_gene_list = list()
subtypes = unique(training_types)
for(cell_type in subtypes){
    marker_gene_list[[cell_type]] = identify_marker_genes(
        training_data,
        training_types,
        cell_type,
        nr_marker_genes = 100
    )
}
# Prepare training data objects
training_mat_bseq = new(
    "ExpressionSet",
    exprs = as.matrix(training_data)
)
fData(training_mat_bseq) = data.frame(rownames(training_data))
pData(training_mat_bseq) = data.frame(training_types)
# Train deconvolution basis matrix
basis = suppressMessages(
    bseqsc_basis(
        training_mat_bseq,
        marker_gene_list,
        clusters = 'training_types',
        samples = colnames(training_data),
        ct.scale = FALSE
    )
)

### Simulate bulk data, process and prepare
alpha_sim_bulk = simulateCellTypes(referenceCellTypes = alpha_bulk,
                                   markerGenes = as.character(unlist(marker_gene_list)),
                                   numSamples = 100)
beta_sim_bulk = simulateCellTypes(referenceCellTypes = beta_bulk,
                                  markerGenes = as.character(unlist(marker_gene_list)),
                                  numSamples = 100)
delta_sim_bulk = simulateCellTypes(referenceCellTypes = delta_bulk,
                                   markerGenes = as.character(unlist(marker_gene_list)),
                                   numSamples = 100)
test_cellTypes = c(rep("Alpha", 100), rep("Beta", 100), rep("Delta", 100))

### Estimate parameters for perturbation from simulated data
### perturbation
interval = seq(0, 1, by = 0.05)
results_perc = c()
results_p = c()
exp_markers = rownames(alpha_sim_bulk) %in% as.character(unlist(marker_gene_list))
exp_markers = cbind(alpha_sim_bulk[exp_markers,], beta_sim_bulk[exp_markers,])
deviation = sd(c(t(exp_markers)))

### Deconvolution of simulated data with increasing noise level
for(j in 1:length(interval)){
    # Alpha expression + noise
    alpha_sim_bulk_noise = c()
    for(i in 1:ncol(alpha_sim_bulk)){
        noise = rnorm(length(alpha_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = alpha_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        alpha_sim_bulk_noise = cbind(alpha_sim_bulk_noise, current)
    }
    # Beta expression + noise
    beta_sim_bulk_noise = c()
    for(i in 1:ncol(beta_sim_bulk)){
        noise = rnorm(length(beta_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = beta_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        beta_sim_bulk_noise = cbind(beta_sim_bulk_noise, current)
    }
    # Delta expression + noise
    delta_sim_bulk_noise = c()
    for(i in 1:ncol(delta_sim_bulk)){
        noise = rnorm(length(delta_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = delta_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        delta_sim_bulk_noise = cbind(delta_sim_bulk_noise, current)
    }
    #combined_sim = cbind(alpha_sim_bulk_noise, beta_sim_bulk_noise)
    combined_sim = cbind(alpha_sim_bulk_noise, beta_sim_bulk_noise, delta_sim_bulk_noise)
    rownames(combined_sim) = rownames(alpha_sim_bulk)
    # Deconvolution
    fit = bseqsc_proportions(
        as.matrix(combined_sim),
        basis, 
        verbose = TRUE,
        log = FALSE,
        perm = 500,
        absolute = TRUE
    )
    # Extract results
    res_coeff = t(fit$coefficients)
    res_coeff[is.na(res_coeff)] = 0.0
    res_cor   = fit$stats
    res_cor[is.na(res_cor)] = 0.0
    # Evaluate performance: p-value <= 0.03, correct cell-type predicted
    predicted_types = apply(res_coeff, 1, evaluate)
    correct = test_cellTypes == predicted_types
    sig_prop = apply(res_coeff, 1, evaluate2)
    correct = correct & sig_prop
    only_one = apply(res_coeff, 1, evaluate3)
    correct = correct & only_one
    p_value = which(res_cor[,1] > 0.03)
    correct[p_value] = FALSE
    # Keep track of mean p-value and % predicted correct
    results_p = c(results_p, mean(res_cor[,1]))
    results_perc = c(results_perc, (sum(correct) / 300))
    
}