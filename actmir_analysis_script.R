# Load PBTA dataset with miRNA and mRNA dataset
met <- readRDS("meta_pbta_miRNA_plus_mRNA.rds")
miRexpr <- readRDS("norm_counts_pbta_mirna.rds")
miRexpr <- as.matrix(miRexpr)
expr <- readRDS("norm_counts_pbta_mrna.rds")
mirna.ids <- readRDS("norm_counts_pbta_mirna_miR-ids.rds")
targetmat <- read.delim('PredictedTarget_dataforBioinformatics.txt',  check.names=F, row.names=1)
target <- rownames(targetmat)[which(!is.na(targetmat[,"MIMAT0018184"]))]
miRNA = c("MIMAT0010195")
cutoff=0.35

# Run actMir analysis
source("./Infer_miRactivity_forBioinformatics.R")
miRact <- InfermiRactivity(miRNA, miRexpr, expr, target, cutoff)

# Export activity scores with PBTA metadata
samples <- met$short_histology
rownames(miRact) <- samples

met2 <- met[order(met$sample_id),] # Order the dataframe based on "sample_id" column
ordered_samples <- miRact[order(names(miRact))] # Order the list based on "sample_id"
met2$activity_score <- ordered_samples

write_csv(met2, "pbta_miRactivity_results.csv")


### EXAMPLE RUN ###
#   source("./Infer_miRactivity_forBioinformatics.R")
# expmat <- read.delim("./data/Expr_dataforBioinformatics.txt", check.names=F, row.names=1)
# miRexpmat <- read.delim("./data/miRexp_dataforBioinformatics.txt", check.names=F, row.names=1)
# targetmat <- read.delim("./data/PredictedTarget_dataforBioinformatics.txt",  check.names=F, row.names=1)
# targets <- rownames(targetmat)[which(!is.na(targetmat[,"MIMAT0000062"]))]
# miRact <- InfermiRactivity("MIMAT0000062", miRexpmat, expmat, targets, 0.35) 





