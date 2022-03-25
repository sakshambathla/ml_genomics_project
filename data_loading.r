library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

TCGAbiolinks:::getProjectSummary("TCGA-LIHC")

query = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")


res = getResults(query)
colnames(res) 
res$sample_type <- as.factor(res$sample_type)
summary(res$sample_type)

filter_query = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = filter_query)

data = GDCprepare(filter_query)

dim(data)
colnames(colData(data))
table(data@colData$vital_status)
table(data@colData$definition)
table(data@colData$tissue_or_organ_of_origin)
table(data@colData$gender)
table(data@colData$race)


dim(assay(data))   
head(assay(data)[,1:5])
saveRDS(object = data,
        file = "data.RDS",
        compress = FALSE)
