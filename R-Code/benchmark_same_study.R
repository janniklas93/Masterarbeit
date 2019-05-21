library(AnnotationDbi)
library(bseqsc)
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/marker_genes.R")

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

# perform benchmark on selected training data
res = benchmark_same_study(study_counts = training_data,
                           study_types = training_types,
                           marker_genes = marker_gene_list)

    
benchmark_same_study = function(study_counts,
                                study_types,
                                marker_genes)
{
    # helper function
    evaluate = function(row){
        ind = which.max(row)
        name = names(row)[ind]
    }
    
    # split data according to cell-type 
    alpha_index = which(study_types %in% "Alpha")
    beta_index = which(study_types %in% "Beta")
    gamma_index = which(study_types %in% "Gamma")
    delta_index = which(study_types %in% "Delta")

    # shuffle data randomly
    alpha_index = alpha_index[sample(length(alpha_index))]
    beta_index = beta_index[sample(length(beta_index))]
    gamma_index = gamma_index[sample(length(gamma_index))]
    delta_index = delta_index[sample(length(delta_index))]

    # create 10 folds for each cell-type
    alpha_folds = cut(seq(1, length(alpha_index)), breaks = 10, labels = FALSE)
    beta_folds = cut(seq(1, length(beta_index)), breaks = 10, labels = FALSE)
    gamma_folds = cut(seq(1, length(gamma_index)), breaks = 10, labels = FALSE)
    delta_folds = cut(seq(1, length(delta_index)), breaks = 10, labels = FALSE)

    alpha_sens = c()
    beta_sens = c()
    gamma_sens = c()
    delta_sens = c()
    
    alpha_ppv = c()
    beta_ppv = c()
    gamma_ppv = c()
    delta_ppv = c()
    
    alpha_spec = c()
    beta_spec = c()
    gamma_spec = c()
    delta_spec = c()

    # do cross-validation
    for(i in 1:10){
        # select current test and training partition for each cell-type
        test_index_alpha  = which(alpha_folds == i)
        train_index_alpha = alpha_index[-test_index_alpha]
        test_index_alpha  = alpha_index[test_index_alpha]
        testData_alpha    = study_counts[, test_index_alpha]
        trainData_alpha   = study_counts[, train_index_alpha]
    
        test_index_beta   = which(beta_folds == i)
        train_index_beta  = beta_index[-test_index_beta]
        test_index_beta   = beta_index[test_index_beta]
        testData_beta     = study_counts[, test_index_beta]
        trainData_beta    = study_counts[, train_index_beta]
    
        test_index_gamma  = which(gamma_folds == i)
        train_index_gamma = gamma_index[-test_index_gamma]
        test_index_gamma  = gamma_index[test_index_gamma]
        testData_gamma    = study_counts[, test_index_gamma]
        trainData_gamma   = study_counts[, train_index_gamma]
    
        test_index_delta  = which(delta_folds == i)
        train_index_delta = delta_index[-test_index_delta]
        test_index_delta  = delta_index[test_index_delta]
        testData_delta    = study_counts[, test_index_delta]
        trainData_delta   = study_counts[, train_index_delta]
    
        # build data sets for training and testing with matching pheno data
        trainSet = cbind(trainData_alpha, trainData_beta, trainData_gamma, trainData_delta)
        testSet = cbind(testData_alpha, testData_beta, testData_gamma, testData_delta)
        train_cellTypes = c(as.character(study_types[train_index_alpha]),
                            as.character(study_types[train_index_beta]),
                            as.character(study_types[train_index_gamma]),
                            as.character(study_types[train_index_delta]))
        train_cellTypes = as.data.frame(train_cellTypes)
        colnames(train_cellTypes) = c("cell_type")
        test_cellTypes = c(as.character(study_types[test_index_alpha]),
                           as.character(study_types[test_index_beta]),
                           as.character(study_types[test_index_gamma]),
                           as.character(study_types[test_index_delta]))
    
        # prepare input for SVR
        single_exp = new("ExpressionSet", exprs = as.matrix(trainSet))
        fData(single_exp) = data.frame(rownames(trainSet))
        pData(single_exp) = train_cellTypes
    
        # Build basis matrix
        B = bseqsc_basis(
            single_exp,
            marker_genes,
            clusters = "cell_type",
            samples = colnames(trainSet),
            ct.scale = FALSE
        )
        fit = bseqsc_proportions(
            as.matrix(testSet),
            B, 
            verbose = TRUE,
            log = FALSE,
            perm = 500,
            absolute = TRUE
        )
    
        # evaluate results
        res_coeff = t(fit$coefficients)
        res_coeff[is.na(res_coeff)] = 0.0
        res_cor   = fit$stats
        res_cor[is.na(res_cor)] = 0.0

        # find out which predictions were not significant and which predicted wrong cell-type
        null = apply(res_coeff, 1, sum)
        null_index = which(null == 0)
        res_coeff = res_coeff[-null_index,]
        test_cellTypes = test_cellTypes[-null_index]
        
        predicted_types = apply(res_coeff, 1, evaluate)
        correct = predicted_types == test_cellTypes
        p_value = which(res_cor[,1] > 0.03)
        correct[p_value] = FALSE
    
        alpha_correct = correct[which(test_cellTypes %in% "Alpha")]
        beta_correct = correct[which(test_cellTypes %in% "Beta")]
        gamma_correct = correct[which(test_cellTypes %in% "Gamma")]
        delta_correct = correct[which(test_cellTypes %in% "Delta")]
        
        not_alpha = which(!test_cellTypes %in% "Alpha")
        not_beta = which(!test_cellTypes %in% "Beta")
        not_gamma = which(!test_cellTypes %in% "Gamma")
        not_delta = which(!test_cellTypes %in% "Delta")
        
        alpha_pred = sum(predicted_types %in% "Alpha")
        beta_pred = sum(predicted_types %in% "Beta")
        gamma_pred = sum(predicted_types %in% "Gamma")
        delta_pred = sum(predicted_types %in% "Delta")
        
        # calculate PPV
        alpha_ppv = c(alpha_ppv, sum(alpha_correct) / alpha_pred)
        beta_ppv = c(beta_ppv, sum(beta_correct) / beta_pred)
        gamma_ppv = c(gamma_ppv, sum(gamma_correct) / gamma_pred)
        delta_ppv = c(delta_ppv, sum(delta_correct) / delta_pred)
        
        # calculate sensitivity
        alpha_sens = c(alpha_sens, sum(alpha_correct) / length(alpha_correct))
        beta_sens = c(beta_sens, sum(beta_correct) / length(beta_correct))
        gamma_sens = c(gamma_sens, sum(gamma_correct) / length(gamma_correct))
        delta_sens = c(delta_sens, sum(delta_correct) / length(delta_correct))
        
        spec_count = c(0, 0, 0, 0)
        names(spec_count) = c("Alpha", "Beta", "Gamma", "Delta")
        # calculate specificity
        for(j in 1:length(test_cellTypes)){
            index = which(colnames(res_coeff) == test_cellTypes[j])
            cell_value = res_coeff[j, index]
            perc = res_coeff[j,] / cell_value
            names_cell = names(perc[perc > 0.5])
            remove = which(names_cell %in% test_cellTypes[j])
            names_cell = names_cell[-remove]
            if(length(names_cell) > 0){
                for(k in 1:length(names_cell)){
                    spec_count[names_cell[k]] = spec_count[names_cell[k]] + 1
                }
            }
            
        }
        alpha_spec = c(alpha_spec, (length(not_alpha) - spec_count["Alpha"]) / length(not_alpha))
        beta_spec = c(beta_spec, (length(not_beta) - spec_count["Beta"]) / length(not_beta))
        gamma_spec = c(gamma_spec, (length(not_gamma) - spec_count["Gamma"]) / length(not_gamma))
        delta_spec = c(delta_spec, (length(not_delta) - spec_count["Delta"]) / length(not_delta))
    }
    
    res = list("alpha_sens" = alpha_sens,
                "beta_sens" = beta_sens,
                "gamma_sens" = gamma_sens,
                "delta_sens" = delta_sens,
                "alpha_ppv" = alpha_ppv,
                "beta_ppv" = beta_ppv,
                "gamma_ppv" = gamma_ppv,
                "delta_ppv" = delta_ppv,
                "alpha_spec" = alpha_spec,
                "beta_spec" = beta_spec,
                "gamma_spec" = gamma_spec,
                "delta_spec" = delta_spec)
    return(res)
}   
