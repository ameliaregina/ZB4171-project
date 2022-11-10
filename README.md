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
GSE200218 - Using RStudio Server Image
GSE181919 - Using RStudio Server Image

## 3. What has already been done?
GSE202695
  - Script Present
  - Generation of plots figures done for population

GSE200218
  - NULL - Progress removed due to server crash

GSE181919
  - NULL - Progress removed due to server crash

EC2 AWS RStudio-Server AMI 
  - Up-to-date RStudio-Server ver., R ver. and Ubuntu

EC2 AWS LouisAslett AMI
  - Working
  - Currently being used to cluster GSE200218 and GSE181919 using Seurat

## 4. What has yet to be done?
GSE202695
  - Data analysis of gene expression and the pathways they regulate
  - Exploration of whole workflow on individual clusters

GSE200218
  - Generation of Plots and figures
  - Data analysis of Gene expression
  - Exploration of pathways regulated by highlighted genes

GSE181919
  - Same as 218

EC2 AWS RStudio-Server AMI 
  - Have not finished installing all necessary packages for workflow
  - Verify functionality of Seurat and Monocle 3

## 5. For each member, please state:
### 1. Other tasks outside of ZB that you need to
### 2. Date of Final Exams
### 3. When you need to start studying for [- 2.]. 
1. Wei Jing
   1-1. Quiz 
        11/11/2022
   1-2. 19/11/2022 and 21/11/2022    
   1-3. Willing to Half-Half until 16/11/2022
   Side-note: if unable to finish during reading week, willing to take over after 21/11/2022

2. Weilin
   2-1. Quiz | Personal Plans | Rushing Report
        11/11/2022 | 12/11/2022 - NOON | 13/11/2022
   2-2. 19/11/2022 and 28/11/2022
   2-3. Willing to Half-Half

3. Christopher
   3-1. Assignment
        12/11/2022
   3-2. 22/11/2022 and 01/12/2022
   3-3. Full-time till 15th - latest 17th Nov

4. Amelia
   4-1. MA Assignment | LSM4242 Report
        12/11/2022 | 13/11/2022
   4-2. 19/11/2022 and 01/12/2022
   4-3. Willing to Half-Half until 16/11/2022

## 6. For items stated in [4.], Please elaborate on:
### 1. Who would be responsible for completing this process?
### 2. What resources are required?
### 3. How long would it take for you to complete the task?
### 4. Provide an estimate of what day you expect the task to be done by?
GSE202695: Proceed with clustering but mention how the dataset might be inappropriate and therefore we 
will refrain from making any assumptions [IC: Chris]

EC2 RStudio AMI: Get all the required dependencies downloaded and ensure all programs work fine by testing 
with whatever dataset is in the S3 bucket. [IC: Wei Jing]

GSE200218 and GSE181919: we would require a suitable EC2 instance with RStudio to run due to the large 
nature of the dataset, however the pipeline will be the same using the original script written
however, the script may require tweaking to suit the formatting of fils in each dataset [IC: Whoever available]

Presentation and Report: Presentation Tuesday, Report Wednesday (BEST CASE SCENARIO) [ICs: Wei Lin, Amelia]


