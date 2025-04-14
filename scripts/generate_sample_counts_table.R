library(SummarizedExperiment)
source("PanClockHelperFunctions.R")

layers <- c("miRNA", "RPPA", "RNAseq", "methylation", "BinaryMutation", "SCNA", "RNAseq_Normal")
sampleCounts <- makeSampleCountsTable(layers)
write.csv(sampleCounts, "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/Eitan/Pancancer_Aging_Clocks/scripts/dataset/allSampleCounts.csv")
