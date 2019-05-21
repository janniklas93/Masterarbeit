library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(EPIC)
library(dplyr)

### Compute immune cell fractions 
meta = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/Meta_information.tsv"
meta = read.table(meta, header = TRUE, sep = "\t")
rownames(meta) = meta$Name

panNEN_data = "~/Documents/Bioinformatik/Master/Masterarbeit/Expressionsdaten/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv"
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
bulk_meta = meta[meta$Name %in% colnames(bulk_data),]
nec_net_vec = as.character(bulk_meta$NEC_NET)
nec_index = which(nec_net_vec %in% "NEC")
nec_samples = rownames(bulk_meta)[nec_index]
net_index = which(nec_net_vec %in% "NET")
net_samples = rownames(bulk_meta)[net_index]

out = EPIC(bulk = bulk_data)
cell_fractions = out$cellFractions
cell_fractions = as.data.frame(cell_fractions)
other = 1 - apply(select(cell_fractions, -c(otherCells)), 1, sum)
cell_fractions$otherCells = other
cell_fractions = cell_fractions * 100

sample_id = data.frame("sample_id" = rownames(cell_fractions))
cell_fractions = cbind(cell_fractions, sample_id)
cell_fractions = melt(cell_fractions, id.vars = "sample_id", measure.vars = c("Bcells", "CAFs", "CD4_Tcells", "CD8_Tcells", "Endothelial", "Macrophages", "NKcells", "otherCells"))
cell_fractions$variable = relevel(cell_fractions$variable, "otherCells")
group = rep("PanNET", nrow(cell_fractions))
group[cell_fractions$sample_id %in% nec_samples] = "PanNEC"
cell_fractions = cbind(cell_fractions, group)
cell_fractions$group = relevel(cell_fractions$group, "PanNET")
