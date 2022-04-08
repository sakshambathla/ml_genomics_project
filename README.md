# ml_genomics_project

**Survival Prediction for Liver Cancer Patients**

Using TCGA gene expression and clinical data, we trained different classifiers for the task of prognosis prediction. We identified different clinical features 
using Kapan-Meir survival analysis and identified deferentially expressed genes using voom and glmfit. After feature extraction, we divided our data into 80:20 
test-train split and compared the performance of various classification algorithms over the data. We were able to show that this outperforms the baseline of selecting 
the majority class. This pipeline can be further extended to different cancer types and based on data availability, different machine learning algorithms can be used.
