
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
library(e1071)
library("MASS")
library("sgd")
library(rpart)
library(class)


limma_res = readRDS(file = "limma_res.RDS")

limma_res$topGenes

plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}
res_pca = plot_PCA(limma_res$voomObj, "definition")

clinical = data
clin_df = clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status", 
                     "days_to_death",
                     "days_to_last_follow_up",
                     "age_at_diagnosis",
                     "race",
                     "gender",
                     "ajcc_pathologic_stage")]


clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,clin_df$days_to_death,clin_df$days_to_last_follow_up)
Surv(clin_df$overall_survival, clin_df$deceased)
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender

fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)
print(fit)
ggsurvplot(fit, data=clin_df, pval = T)
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, risk.table.col="strata")



fit2 = survfit(Surv(overall_survival, deceased) ~ ajcc_pathologic_stage, data=clin_df)
print(fit2)
ggsurvplot(fit2[c(1,2,4)], data=clin_df, pval=T, risk.table=T, risk.table.col="strata")


fit3 = survfit(Surv(overall_survival, deceased) ~ race, data=clin_df)
print(fit3)
ggsurvplot(fit3[c(2,3,5)], data=clin_df, pval=T, risk.table=T, risk.table.col="strata")


d_mat = as.matrix(t(limma_res$voomObj$E))

d_resp = as.factor(limma_res$voomObj$targets$definition)

set.seed(50)

train_ids = createDataPartition(d_resp, p=0.75, list=FALSE)

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

y_train = d_resp[train_ids]
y_test  = d_resp[-train_ids]

res = cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = 0.5,
  family = "binomial"
)

y_pred = predict(res, newx=x_test, type="class", s="lambda.min")
confusion_matrix = table(y_pred, y_test)

print(confusion_matrix)

print(paste0("Sensitivity: ",sensitivity(confusion_matrix)))

print(paste0("Specificity: ",specificity(confusion_matrix)))

print(paste0("Precision: ",precision(confusion_matrix)))


res_coef = coef(res, s="lambda.min") 

res_coef = res_coef[res_coef[,1] != 0,]
length(res_coef)

res_coef

relevant_genes = names(res_coef) 
relevant_genes = relevant_genes[-1]

head(relevant_genes) # few select genes
head(limma_res$voomObj$genes)
relevant_gene_names = limma_res$voomObj$genes[relevant_genes,"external_gene_name"]

head(relevant_gene_names) # few select genes (with readable names now)


print(intersect(limma_res$topGenes$ensembl_gene_id, relevant_genes))

hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)


clust = function(x) hclust(x, method="complete")

dist = function(x) as.dist((1-cor(t(x)))/2)


colorLimmaGenes = ifelse(
  
  (relevant_genes %in% limma_res$topGenes$ensembl_gene_id),
  "green", 
  "white" 
)

gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          
  density.info="none",  
  trace="none",         
  col=hmcol,            
  labRow=relevant_gene_names, 
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         
  ColSideColors=as.character(as.numeric(d_resp)), 
  dendrogram="both",    
  hclust = clust,       
  distfun = dist,       
  cexRow=.6,            
  margins=c(1,5)        
)



expr_df = limma_res$topGenes

# print the first row, to see the gene name, the logFC value and the p-value
print(expr_df[1:88, ])

# get the ensembl gene id of the first row
gene_ids = expr_df[1:88, "ensembl_gene_id"]

# also get the common gene name of the first row
gene_names = expr_df[1:88, "external_gene_name"]

# we now have selected a gene.
# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)

expr_diseased = d_mat[rownames(clin_df), gene_ids[1]]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_ids[1]]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased gene 1", "Healthy gene 1"), main="Distribution of gene expression")

expr_diseased = d_mat[rownames(clin_df), gene_ids[2]]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_ids[2]]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased gene 2", "Healthy gene 2"), main="Distribution of gene expression")

expr_diseased = d_mat[rownames(clin_df), gene_ids[3]]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_ids[3]]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased gene 3", "Healthy gene 3"), main="Distribution of gene expression")

# get the expression values for the top de genes
i = 1
for(gene_name in gene_names) {
  clin_df[gene_name] = d_mat[rownames(clin_df), gene_ids[i]]
  i = i+1
}


clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)


clin_df$long_term = ifelse(clin_df$overall_survival >= 365, 1, 0)

colnames(clin_df) 

clin_df <- subset (clin_df, select = -c(patient,vital_status, days_to_death, days_to_last_follow_up,overall_survival, deceased))




encode_ordinal <- function(x, order = unique(x)) {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x
}


clin_df$gender = encode_ordinal(clin_df$gender)
clin_df$race = encode_ordinal(clin_df$race)
clin_df$ajcc_pathologic_stage = encode_ordinal(clin_df$ajcc_pathologic_stage)
clin_df <- na.omit(clin_df)
scaled_df = clin_df
scaled_df['age_at_diagnosis'] = scale(scaled_df['age_at_diagnosis'])

split1<- sample(c(rep(0, 0.8 * nrow(scaled_df)), rep(1, 0.2 * nrow(scaled_df))))
train <- scaled_df[split1 == 0, ]            
test <- scaled_df[split1== 1, ]    

precision <- function(matrix) {
  # True positive
  tp <- matrix[2, 2]
  # false positive
  fp <- matrix[1, 2]
  return (tp / (tp + fp))
}

recall <- function(matrix) {
  # true positive
  tp <- matrix[2, 2]# false positive
  fn <- matrix[2, 1]
  return (tp / (tp + fn))
}



##Logistic Regression
mylogit = glm(long_term ~ ., data = train, family = "binomial")
predict <- predict(mylogit, test, type = 'response', method = 'looCV')
table_mat <- table(test$long_term, predict > 0.5)
accuracy_test <- sum(diag(table_mat)) / sum(table_mat)
prec <- precision(table_mat)
rec <- recall(table_mat)
f1 <- 2 * ((prec * rec) / (prec + rec))
## LR
table_mat
accuracy_test
prec
rec
f1


## SVM
classifier = svm(formula = long_term ~ .,data = train,type = 'C-classification',kernel = 'linear')
y_pred = predict(classifier, newdata = test, method = 'looCV')
table_mat = table(test$long_term, y_pred)

accuracy_test <- sum(diag(table_mat)) / sum(table_mat)
prec <- precision(table_mat)
rec <- recall(table_mat)
f1 <- 2 * ((prec * rec) / (prec + rec))
## SVM
table_mat
accuracy_test
prec
rec
f1

## LDA
ldamodel = lda(long_term ~ .,data = train)
y_pred = predict(ldamodel, newdata = test)
table_mat = table(test$long_term, y_pred['class'][1]$class)
accuracy_test <- sum(diag(table_mat)) / sum(table_mat)
prec <- precision(table_mat)
rec <- recall(table_mat)
f1 <- 2 * ((prec * rec) / (prec + rec))
## LDA
table_mat
accuracy_test
prec
rec
f1

## Decision Trees
dtr <- rpart(long_term~., data =train, method = 'class')
y_pred = predict(dtr, newdata = test , method = 'looCV')
ans = y_pred[, 2] - y_pred[, 1] 
ans = ifelse(ans>0,1,0)
table_mat = table(test$long_term, ans)
accuracy_test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_test
prec <- precision(table_mat)
rec <- recall(table_mat)
f1 <- 2 * ((prec * rec) / (prec + rec))
## decision tree
table_mat
accuracy_test
prec
rec
f1

## KNN
pr <- knn(train,test,cl=train$long_term,k=13)
table_mat <- table(pr,test$long_term)
accuracy_test <- sum(diag(table_mat)) / sum(table_mat)
prec <- precision(table_mat)
rec <- recall(table_mat)
f1 <- 2 * ((prec * rec) / (prec + rec))
## KNN
table_mat
accuracy_test
prec
rec
f1
