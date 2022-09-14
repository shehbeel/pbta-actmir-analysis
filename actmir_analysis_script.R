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

mirnas <- row.names(miRexpr)
df <- data.frame(matrix(nrow = 244, ncol=0))
row.names(df) <- colnames(miRexpr)
counter <- 0 # Using a counter to let me know when each loop is complete, so I can have runtime reference

for (i in mirnas)
{
  miRact <- InfermiRactivity(i, miRexpr, expr, target, cutoff)
  #print(miRact)
  df[i] <- cbind(miRact, df)
  counter = counter + 1
  print(counter)
}

final.df <- cbind(met, df)

# Export inferred miRNA Activity Scores dataframe
write.csv(final.df, "pbta_miRactivity_results.csv", row.names=FALSE)



### EXAMPLE RUN ###
#   source("./Infer_miRactivity_forBioinformatics.R")
# expmat <- read.delim("./data/Expr_dataforBioinformatics.txt", check.names=F, row.names=1)
# miRexpmat <- read.delim("./data/miRexp_dataforBioinformatics.txt", check.names=F, row.names=1)
# targetmat <- read.delim("./data/PredictedTarget_dataforBioinformatics.txt",  check.names=F, row.names=1)
# targets <- rownames(targetmat)[which(!is.na(targetmat[,"MIMAT0000062"]))]
# miRact <- InfermiRactivity("MIMAT0000062", miRexpmat, expmat, targets, 0.35) 





