library(bseqsc)
library(ggplot2)
library(stringr)
library(reshape2)
library(RColorBrewer)
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/marker_genes.R")

single_cell_seq = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Merged_Expression_Data/Baron_Diff_Haber_Hisc_merged.tsv"
panNEN_data = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv"
meta_data = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Meta_information.tsv"
training_nr_marker_genes = 100
nr_permutations = 500

deconvolution_panNEN = function(
    single_cell_seq = "",
    panNEN_data = "",
    meta_data = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Meta_information.tsv",
    training_nr_marker_genes = 100,
    nr_permutations = 100
){

    ### Read in data
    meta = read.table(meta_data, header = TRUE, sep = "\t")
    rownames(meta) = meta$Name
    training_data = read.table(single_cell_seq,
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE
    )
    colnames(training_data) = str_replace(colnames(training_data), pattern = "\\.", "_")
    colnames(training_data) = str_replace_all(colnames(training_data), pattern = "^X", "")
    training_types = as.character(meta$Subtype[rownames(meta) %in% colnames(training_data)])
    
    # Read in PanNEN bulk data
    bulk_data = read.table(
                    panNEN_data,
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)
    colnames(bulk_data) = stringr::str_replace_all(
        colnames(bulk_data),
        pattern = "^X",
        ""
    )
    # Process meta data and extract PanNEN sample classification
    bulk_meta = meta[meta$Name %in% colnames(bulk_data),]
    nec_net_vec = as.character(bulk_meta$NEC_NET)
    nec_index = which(nec_net_vec %in% "NEC")
    nec_samples = rownames(bulk_meta)[nec_index]
    net_index = which(nec_net_vec %in% "NET")
    net_samples = rownames(bulk_meta)[net_index]

    ### Data cleaning
    row_var = apply(training_data, FUN = var, MARGIN = 1)
    training_data = training_data[row_var != 0,]
    training_data = training_data[rowSums(training_data) >= 1,]

    ### Determine marker gene by differential expression
    marker_gene_list = list()
    subtypes = c("Alpha", "Beta", "Gamma", "Delta", "HISC")
    for(cell_type in subtypes){
        marker_gene_list[[cell_type]] = identify_marker_genes(
            training_data,
            training_types,
            cell_type,
            training_nr_marker_genes
        )
    }

    ### Prepare training
    training_mat_bseq = new(
        "ExpressionSet",
        exprs = as.matrix(training_data)
    )
    fData(training_mat_bseq) = data.frame(rownames(training_data))
    pData(training_mat_bseq) = data.frame(training_types)

    basis = suppressMessages(
        bseqsc_basis(
            training_mat_bseq,
            marker_gene_list,
            clusters = 'training_types',
            samples = colnames(training_data),
            ct.scale = FALSE
        )
    )
    
    ### Deconvolution with trained model
    deconvolution_data = new("ExpressionSet", exprs=as.matrix(bulk_data))
    fit = bseqsc_proportions(
            deconvolution_data,
            basis,
            verbose = FALSE,
            absolute = TRUE,
            log = FALSE,
            perm = nr_permutations,
            QN = FALSE
        )

    res_coeff = t(fit$coefficients)
    res_coeff[is.na(res_coeff)] = 0.0
    res_cor   = fit$stats
    res_cor[is.na(res_cor)] = 0.0
    
    colnames(res_coeff) = c("Alpha", "Beta", "Gamma", "Delta", "Progenitor", "HISC")
    training_types = c("Alpha", "Beta", "Gamma", "Delta", "Progenitor", "HISC")
    
    results = cbind(res_coeff, res_cor)
    return(results)
}