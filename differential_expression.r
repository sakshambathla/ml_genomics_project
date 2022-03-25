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



data = readRDS(file = "data.RDS")

pipeline = function(
    data,
    condition_variable,
    reference_group=NULL){
  
  design_factor = colData(data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(data),
                samples=colData(data),
                genes=as.data.frame(rowData(data)))
  
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  
  return(
    list(
      voomObj=v, 
      fit=fit, 
      topGenes=topGenes
    )
  )
}

limma_res = pipeline(
  data=data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)

saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)
