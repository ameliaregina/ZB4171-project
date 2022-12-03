# ZB4171: Characterizing Phenotypic Plasticity in Metastatic Cell Trajectories (Project Timeline)

### **Group members: Amelia, Chris, Wei Jing, Wei Lin**

## 1. What our goal is?
### - Please provide details on the workflow and expected results or figures or tables
To conduct preliminary studies on pan-cancer primary and metastatic tumors to characterize phenotypic plasticity by cell trajectory inferences and differential gene expression analyses.
By looking into various cancer datasets we hope to identify potential genes that may be responsible for initiating the hallmark of phenotypic plasticity in cancer. Such that we can suggest novel metastasis biomarkers or therapeutic targets for prevention.

### - Current issue: 
GSE202695 (breast primary + lung metastasis) was our proposed dataset, but upon generating the results as we initially intended, we realized that the way the original authors generated the dataset may not be appropriate for cell trajectory inference; hence we chose to analyze 218 and 1919 as the abstract suggests data collected
is appropriate for this workflow.

## 2. What datasets are we using?
GSE202695 - Done at home computer
GSE181919 - Using RStudio Server on EC2

## 3. What was done?
GSE202695
  - Whole dataset clustering + trajectory (+3D trajectory)
  - Per model type clustering trajectory 

GSE181919
  - Whole dataste clustering + trajectory
  - Clustering + trajectory of just CA (primary cancer) and LN (lymph node metastasis)
