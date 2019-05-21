library(AnnotationDbi)
library(bseqsc)
library(stringr)
source("~/Documents/Bioinformatik/Master/Masterarbeit/R_Scripts/marker_genes.R")

meta_data = read.table("~/Downloads/ArtDeco/inst/Data/Meta_information/Meta_information.tsv", header = TRUE, sep = "\t")

# read in HISC (Haber)
haber_counts = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Human_Mouse_HSC/Haber.tsv"
haber_counts = read.table(haber_counts , sep ="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colnames(haber_counts) = str_replace_all(colnames(haber_counts) , pattern = "^X", "")
rownames(haber_counts) = str_to_upper(rownames(haber_counts))
rownames(haber_counts)[rownames(haber_counts) == "INS-IGF2" ] = "INS"
haber_meta = meta_data[meta_data$Sample_ID %in% colnames(haber_counts),]
haber_counts = haber_counts[, haber_meta$Subtype %in% c("HISC")]
haber_types = rep("hisc", ncol(haber_counts))

# read in Prog (Stanescu)
stanescu_counts = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Mouse_progenitor_pancreas_scRNA/Stanescu.tsv"
stanescu_counts = read.table(stanescu_counts , sep = "\t", header = T, stringsAsFactors = F)
rownames(stanescu_counts) = str_to_upper(rownames(stanescu_counts))
rownames(stanescu_counts)[rownames(stanescu_counts) == "INS1"] = "INS"
stanescu_types = rep("prog", ncol(stanescu_counts))

### integrate
merge_genes = intersect(rownames(haber_counts),rownames(stanescu_counts))
new_mat = as.data.frame(
    cbind(
        haber_counts[rownames(haber_counts) %in% merge_genes,],
        stanescu_counts[rownames(stanescu_counts) %in% merge_genes,]
    )
)
rownames(new_mat) = merge_genes
cell_types = c(haber_types, stanescu_types)

prog_new = identify_marker_genes(new_mat, as.character(cell_types), subtype = "prog")
hisc_new = identify_marker_genes(new_mat, as.character(cell_types), subtype = "hisc")
pancreas_markers = list(
    "hisc" = hisc_new,
    "prog" = prog_new
)

# helper function
evaluate = function(row){
    ind = which.max(row)
    name = names(row)[ind]
}

# split data according to cell-type 
hisc_index = which(cell_types %in% "hisc")
prog_index = which(cell_types %in% "prog")

# shuffle data randomly
hisc_index = hisc_index[sample(length(hisc_index))]
prog_index = prog_index[sample(length(prog_index))]

# create 10 folds for each cell-type
hisc_folds = cut(seq(1, length(hisc_index)), breaks = 10, labels = FALSE)
prog_folds = cut(seq(1, length(prog_index)), breaks = 10, labels = FALSE)

pancreas_markers$hisc = pancreas_markers$hisc[(
    pancreas_markers$hisc %in% rownames(new_mat))]
pancreas_markers$prog = pancreas_markers$prog[(
    pancreas_markers$prog %in% rownames(new_mat))]

hisc_sens = c()
hisc_ppv = c()
hisc_spec = c()

prog_sens = c()
prog_ppv = c()
prog_spec = c()

# do cross-validation
for(i in 1:10){
    # select current test and training partition for each cell-type
    test_index_hisc  = which(hisc_folds == i)
    train_index_hisc = hisc_index[-test_index_hisc]
    test_index_hisc  = hisc_index[test_index_hisc]
    testData_hisc    = new_mat[, test_index_hisc]
    trainData_hisc   = new_mat[, train_index_hisc]
    
    test_index_prog = which(prog_folds == i)
    train_index_prog = prog_index[-test_index_prog]
    test_index_prog = prog_index[test_index_prog]
    testData_prog = new_mat[, test_index_prog]
    trainData_prog = new_mat[, train_index_prog]
    
    # build data sets for training and testing with matching pheno data
    trainSet = cbind(trainData_hisc, trainData_prog)
    testSet = cbind(testData_hisc, testData_prog)
    train_cellTypes = c(as.character(cell_types[train_index_hisc]),
                        as.character(cell_types[train_index_prog]))
    train_cellTypes = as.data.frame(train_cellTypes)
    colnames(train_cellTypes) = c("cell_type")
    test_cellTypes = c(as.character(cell_types[test_index_hisc]),
                       as.character(cell_types[test_index_prog]))
    
    # prepare input for CIBERSORT
    single_exp = new("ExpressionSet", exprs = as.matrix(trainSet))
    fData(single_exp) = data.frame(rownames(trainSet))
    pData(single_exp) = train_cellTypes
    
    # Run CIBERSORT
    B = bseqsc_basis(
        single_exp,
        pancreas_markers,
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
    res_coeff[res_cor[,"Correlation"] <= .3, ] = .0
    res_cor[ res_cor[,"Correlation"] <= .3, "Correlation" ] = 0.0
    
    # find out which predictions were not significant and which predicted wrong cell-type
    predicted_types = apply(res_coeff, 1, evaluate)
    correct = predicted_types == test_cellTypes
    p_value = which(res_cor[,1] > 0.03)
    correct[p_value] = FALSE
    
    hisc_correct = correct[which(test_cellTypes %in% "hisc")]
    prog_correct = correct[which(test_cellTypes %in% "prog")]
    
    not_hisc = which(!test_cellTypes %in% "hisc")
    not_prog = which(!test_cellTypes %in% "prog")
    
    hisc_pred = sum(predicted_types %in% "hisc")
    prog_pred = sum(predicted_types %in% "prog")
    
    # calculate PPV
    hisc_ppv = c(hisc_ppv, sum(hisc_correct) / hisc_pred)
    prog_ppv = c(prog_ppv, sum(prog_correct) / prog_pred)
    
    # calculate sensitivity
    hisc_sens = c(hisc_sens, sum(hisc_correct) / length(hisc_correct))
    prog_sens = c(prog_sens, sum(prog_correct) / length(prog_correct))
    
    # calculate specificty
    hisc_spec = c(hisc_spec, sum(!predicted_types[not_hisc] %in% "hisc") / length(not_hisc))
    prog_spec = c(prog_spec, sum(!predicted_types[not_prog] %in% "prog") / length(not_prog))
}

res = list("hisc_sens" = hisc_sens,
           "hisc_ppv"  = hisc_ppv,
           "hisc_spec" = hisc_spec,
           "prog_sens" = prog_sens,
           "prog_ppv"  = prog_ppv,
           "prog_spec" = prog_spec)

