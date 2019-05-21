library(AnnotationDbi)
library(bseqsc)
library(ggplot2)
library(OptimalCutpoints)
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/marker_genes.R")
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/Simulation.R")

alpha_beta = read.table(
    "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Purified/alpha_beta_purified/purified_alpha_beta_bulk_tpm.tsv",
    header = TRUE,
    sep = "\t"
)
alpha_bulk = alpha_beta[,grep("Alpha", colnames(alpha_beta))]
beta_bulk = alpha_beta[,grep("Beta", colnames(alpha_beta))]

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
training_data = lawlor_counts
training_types = lawlor_types
colnames(training_data) = str_replace(colnames(training_data), pattern = "\\.", "_")
colnames(training_data) = str_replace_all(colnames(training_data), pattern = "^X", "")
# Data cleaning
row_var = apply(training_data, FUN = var, MARGIN = 1)
training_data = training_data[row_var != 0,]
training_data = training_data[rowSums(training_data) >= 1,]
# Determine marker genes by differential expression
marker_gene_list = list()
subtypes = c("Alpha", "Beta", "Gamma", "Delta")
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
negative_controls = simulateNegativeControls(nGenes = nrow(alpha_bulk),
                                             numSamples = 100)
rownames(negative_controls) = rownames(alpha_bulk)

status = c(rep(1, 200), rep(0, 100))
label = c(rep("positive", 200), rep("negative", 100))

### cutoff p-value
interval = c(0, 0.25, 0.5, 0.75, 1)
results_cut = c()
results_sens = c()
results_spec = c()
exp_markers = rownames(alpha_sim_bulk) %in% as.character(unlist(marker_gene_list))
exp_markers = cbind(alpha_bulk[exp_markers,], beta_bulk[exp_markers,])
deviation = sd(c(t(exp_markers)))

for(j in 1:length(interval)){
    # alpha
    alpha_sim_bulk_noise = c()
    for(i in 1:ncol(alpha_sim_bulk)){
        noise = rnorm(length(alpha_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = alpha_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        alpha_sim_bulk_noise = cbind(alpha_sim_bulk_noise, current)
    }
    # beta
    beta_sim_bulk_noise = c()
    for(i in 1:ncol(beta_sim_bulk)){
        noise = rnorm(length(beta_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = beta_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        beta_sim_bulk_noise = cbind(beta_sim_bulk_noise, current)
    }
    negative_bulk_noise = c()
    for(i in 1:ncol(negative_controls)){
        noise = rnorm(length(negative_controls[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = negative_controls[,i] + noise
        index = which(current < 0)
        current[index] = 0
        negative_bulk_noise = cbind(negative_bulk_noise, current)
    }
    combined_sim = cbind(alpha_sim_bulk_noise, beta_sim_bulk_noise, negative_bulk_noise)
    rownames(combined_sim) = rownames(alpha_sim_bulk)

    fit = bseqsc_proportions(
        as.matrix(combined_sim),
        basis, 
        verbose = TRUE,
        log = FALSE,
        perm = 500,
        absolute = TRUE
    )
    
    res_coeff = t(fit$coefficients)
    res_coeff[is.na(res_coeff)] = 0.0
    res_cor   = fit$stats
    res_cor[is.na(res_cor)] = 0.0
    
    prediction = data.frame("Status" = status, "Label" = label, "P_value" = res_cor[,1])
    
    cutoff = optimal.cutpoints(X = "P_value",
                               status = "Status",
                               tag.healthy = 0,
                               methods = "Youden",
                               data = prediction,
                               conf.level = 0.95)
    optimal_cut = cutoff$Youden$Global$optimal.cutoff$cutoff
    results_cut = c(results_cut, optimal_cut)
    optimal_sens = cutoff$Youden$Global$optimal.cutoff$Se
    results_sens = c(results_sens, optimal_sens)
    optimal_spec = cutoff$Youden$Global$optimal.cutoff$Sp
    results_spec = c(results_spec, optimal_spec)
    
}

### cutoff proportions
interval = c(0, 0.25, 0.5, 0.75, 1)
results_cut_prop = c()
results_sens_prop = c()
results_spec_prop = c()
exp_markers = rownames(alpha_sim_bulk) %in% as.character(unlist(marker_gene_list))
exp_markers = cbind(alpha_bulk[exp_markers,], beta_bulk[exp_markers,])
deviation = sd(c(t(exp_markers)))

for(j in 1:length(interval)){
    # alpha
    alpha_sim_bulk_noise = c()
    for(i in 1:ncol(alpha_sim_bulk)){
        noise = rnorm(length(alpha_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = alpha_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        alpha_sim_bulk_noise = cbind(alpha_sim_bulk_noise, current)
    }
    # beta
    beta_sim_bulk_noise = c()
    for(i in 1:ncol(beta_sim_bulk)){
        noise = rnorm(length(beta_sim_bulk[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = beta_sim_bulk[,i] + noise
        index = which(current < 0)
        current[index] = 0
        beta_sim_bulk_noise = cbind(beta_sim_bulk_noise, current)
    }
    negative_bulk_noise = c()
    for(i in 1:ncol(negative_controls)){
        noise = rnorm(length(negative_controls[,i]), 0, interval[j] * log2(deviation))
        noise = 2 ^ noise
        current = negative_controls[,i] + noise
        index = which(current < 0)
        current[index] = 0
        negative_bulk_noise = cbind(negative_bulk_noise, current)
    }
    combined_sim = cbind(alpha_sim_bulk_noise, beta_sim_bulk_noise, negative_bulk_noise)
    rownames(combined_sim) = rownames(alpha_sim_bulk)
    
    fit = bseqsc_proportions(
        as.matrix(combined_sim),
        basis, 
        verbose = TRUE,
        log = FALSE,
        perm = 500,
        absolute = TRUE
    )
    
    res_coeff = t(fit$coefficients)
    res_coeff[is.na(res_coeff)] = 0.0
    res_cor   = fit$stats
    res_cor[is.na(res_cor)] = 0.0
    
    proportions_alpha = res_coeff[1:100,1]
    proportions_beta = res_coeff[101:200,2]
    proportions_neg = apply(res_coeff[201:300,], 1, max)
    
    proportions_combined = c(proportions_alpha, proportions_beta, proportions_neg)
    
    prediction = data.frame("Status" = status, "Label" = label, "Prop" = proportions_combined)
    
    cutoff = optimal.cutpoints(X = "Prop",
                               status = "Status",
                               tag.healthy = 0,
                               methods = "Youden",
                               data = prediction,
                               conf.level = 0.95)
    optimal_cut = cutoff$Youden$Global$optimal.cutoff$cutoff
    results_cut_prop = c(results_cut_prop, optimal_cut)
    optimal_sens = cutoff$Youden$Global$optimal.cutoff$Se
    results_sens_prop = c(results_sens_prop, optimal_sens)
    optimal_spec = cutoff$Youden$Global$optimal.cutoff$Sp
    results_spec_prop = c(results_spec_prop, optimal_spec)
    
}
