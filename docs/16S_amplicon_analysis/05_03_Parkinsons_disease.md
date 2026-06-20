# **Practice 2: Parkinson's disease**
>This chapter contains a training data for 16S amplicon analysis and an example of analysis.

## **Introduction**

In this study, the gut microbiome composition in individuals with Parkinson’s Disease and Healthy Controls was investigated. The goal was to identify taxonomic shifts and associations with disease. High-throughput sequencing data and various bioinformatics tools were employed to analyze microbial diversity and abundance.<br>

Data obtained from the article [«_Parkinson’s Disease and PD Medications Have Distinct Signatures of the Gut Microbiome_»](https://pmc.ncbi.nlm.nih.gov/articles/PMC5469442/) was used. These data include information on microbial samples from healthy individuals and patients with Parkinson’s's disease.<br>

## **Instruction**

You can run commands below in your `RStudio` in `R script`.<br>
Or if you want to write a beautiful & convenient to read laboratory journal you can use `R Markdown`.<br>

----------------------------------------------

## **Step 1: Data Loading and Preprocessing**

Install or call libraries<br>

```r
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, ROCR, PERMANOVA, plyr, data.table, ggplot2, e1071, randomForest, caret, NearestBalance, zCompositions, selbal)
```

???+ note
    If you experience any troubles in installing `NearestBalance` & `selbal` packages please go to [previous chapter](../05_02_Crohns_disease)
 
Then set the working directory.<br>
 
**_Input_**

```r
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)
```

Download the data to work with.<br>

**_Input_**

```r
url <- "https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/data/05_16S_amplicon_analysis/05_03_Parkinsons_disease.zip"

zipF<- "05_03_Parkinsons_disease.zip"

download.file(url, zipF)

outDir<-"."

unzip(zipF,exdir=outDir)

if (file.exists(zipF)) {
  file.remove(zipF)
}
```

Load data from various files.<br>
- `Parkinson_otu_table_L6.tsv` is an OTU (Operational Taxonomic Units) table at level 6 taxonomy (genus level).<br>
- `Parkinson_chao1.tsv` contains alpha diversity metrics (Chao1 richness estimate).<br>
- `sample_info.txt` is metadata, which contains sample information, such as case/control labels.<br>

**_Input_**

```r
biomeData <- read_delim('data/Parkinson_otu_table_L6.tsv',  "\t", col_names = TRUE, escape_double = FALSE, trim_ws = TRUE)
alphaDivData <- read_delim('data/Parkinson_chao1.tsv',  "\t", col_names = TRUE,escape_double = FALSE, trim_ws = TRUE)
metaData <- read_delim('data/sample_info.txt',  "\t", col_names = TRUE, escape_double = FALSE, trim_ws = TRUE)
```

Check if data is loaded correctly as a data frame<br>

**_Input_**

```r
is.data.frame(biomeData)
```

**_Output_**

```
[1] TRUE
```

Ensure that the first column name (`sample_name`) in all datasets matches for merging later<br>

**_Input_**

```r
colnames(biomeData)[1]<-colnames(metaData)[1]
colnames(alphaDivData)[1]<-colnames(metaData)[1]
```

Find the common sample names between the three datasets (`biomeData`, `alphaDivData`, and `metaData`)<br>

**_Input_**

```r
commonSamples<-Reduce(intersect, list(biomeData$sample_name,alphaDivData$sample_name,metaData$sample_name))
```

Set the row names of the data frames to the sample names for easier subsetting later<br>

**_Input_**

```r
rownames(biomeData)<-biomeData$sample_name
rownames(alphaDivData)<-alphaDivData$sample_name
rownames(metaData)<-metaData$sample_name
```

Subset the data to include only the common samples across the datasets<br>

**_Input_**

```r
alphaDivDataS<-alphaDivData[commonSamples,]
biomeDataS<-biomeData[commonSamples,]
metaDataS<-metaData[commonSamples,]
```

Remove the first column (sample name) from biomeDataS as it's redundant now<br>

**_Input_**

```r
biomeDataS<-biomeDataS[,-1]
```

Reassign row names based on sample names from metadata (just to ensure they match)<br>

**_Input_**

```r
rownames(biomeDataS)<-metaDataS$sample_name
```

----------------------------------------------

## **Step 2: Statistical Tests (Normality, U-Test, and GLM)**

### **Normality test**

Shapiro-Wilk test to check if the Chao1 alpha diversity values follow a normal distribution<br>

**_Input_**

```r
shapiro.test(alphaDivDataS$mean_chao1)
```

**_Output_**

```
	Shapiro-Wilk normality test

data:  alphaDivDataS$mean_chao1
W = 0.98753, p-value = 0.0146
```

The data is not normally distributed, the p-value is very small<br>

### **Wilcoxon rank-sum test**

Wilcoxon rank-sum test (non-parametric test) to compare cases and controls<br>
Extract sample names for controls and cases based on `case_control` metadata<br>

**_Input_**

```r
controls<-metaDataS[which(metaDataS$case_control == 'Control'),'sample_name']
cases<-metaDataS[which(metaDataS$case_control == 'Case'),'sample_name']
cases<-cases[[1]]
controls<-controls[[1]]
```

Initialize a matrix to store Wilcoxon test results for each taxonomic group<br>

**_Input_**

```r
wilcoxRes<-matrix(ncol=3,nrow=0)
```

Loop through each taxonomic group (each column of biomeDataS) and perform Wilcoxon test between cases and controls<br>

**_Input_**

```r
for (i in colnames(biomeDataS))
{
  wt<-wilcox.test(biomeDataS[cases,i][[1]],biomeDataS[controls,i][[1]])
  wilcoxRes<-rbind(wilcoxRes, c(i,wt$statistic, wt$p.value))
}
```

Adjust p-values using the False Discovery Rate (FDR) method to control for multiple comparisons<br>

**_Input_**

```r
wicoxPvalAdj<-p.adjust(wilcoxRes[,3], method = 'fdr')
```

Combine the Wilcoxon test results with adjusted p-values into a data frame<br>

**_Input_**

```r
wilcoxRes<-cbind(wilcoxRes, wicoxPvalAdj)
colnames(wilcoxRes)<-c('tax','stat','pval','pval_adj')
wilcoxRes<-as.data.frame(wilcoxRes)
```

Convert p-values to numeric and round them to 3 decimal places for readability<br>

**_Input_**

```r
wilcoxRes$pval<- round(as.numeric(as.character(wilcoxRes$pval)), 3)
wilcoxRes$pval_adj<- round(as.numeric(as.character(wilcoxRes$pval_adj)), 3)
```

Extract significant results where the adjusted p-value is less than 0.05<br>

**_Input_**

```r
wilcoxResSign<-wilcoxRes[which(wilcoxRes$pval_adj <0.05),]
wilcoxResSign
```

**_Output_**

```
tax
4                  Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Mobiluncus
18        Bacteria;Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium
23            Bacteria;Actinobacteria;Actinobacteria;Coriobacteriales;Coriobacteriaceae;Gordonibacter
33                  Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Porphyromonas
64         Bacteria;Firmicutes;Clostridia;Clostridiales;(Eubacteriaceae/Lachnospiraceae);unclassified
65       Bacteria;Firmicutes;Clostridia;Clostridiales;(Eubacteriaceae/Lachnospiraceae)_2;unclassified
75  Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales_Family_XI._Incertae_Sedis;Anaerococcus
81                            Bacteria;Firmicutes;Clostridia;Clostridiales;Eubacteriaceae;Eubacterium
84      Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;(Lachnoclostridium/unclassified)
85                          Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anaerostipes
96                             Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia
211      Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
     stat  pval pval_adj
4   11589 0.000    0.003
18  12781 0.000    0.002
23  11125 0.002    0.042
33  12098 0.000    0.009
64   7651 0.002    0.042
65   7098 0.000    0.004
75  11803 0.001    0.040
81   7436 0.001    0.020
84   7313 0.000    0.009
85   6923 0.000    0.002
96   6692 0.000    0.002
211 12039 0.001    0.035
```

Perform Wilcoxon test for alpha diversity (Chao1) between cases and controls<br>

**_Input_**

```r
wt<-wilcox.test(alphaDivDataS[which(alphaDivDataS$sample_name %in% cases),'mean_chao1'][[1]], alphaDivDataS[which(alphaDivDataS$sample_name %in% controls),'mean_chao1'][[1]])
wt
```

**_Output_**

```
	Wilcoxon rank sum test with continuity correction

data:  alphaDivDataS[which(alphaDivDataS$sample_name %in% cases), "mean_chao1"][[1]] and alphaDivDataS[which(alphaDivDataS$sample_name %in% controls), "mean_chao1"][[1]]
W = 10352, p-value = 0.4215
alternative hypothesis: true location shift is not equal to 0
```

### **Generalized Linear Model (GLM)**

GLM to test associations between microbiome features and metadata (case/control, BMI, sex, age)<br>

Join the metadata with biome data<br>

**_Input_**

```r
glmDF<-inner_join(metaDataS[,c('sample_name','case_control','sex','age','bmi')],biomeData, by = 'sample_name')
```

Filter out samples with missing BMI, age, or sex data<br>

**_Input_**

```r
glmDF<-glmDF[which(!is.na(glmDF$bmi)),]
glmDF<-glmDF[which(!is.na(glmDF$age)),]
glmDF<-glmDF[which(!is.na(glmDF$sex)),]
```

Convert age and BMI columns to numeric type<br>

**_Input_**

```r
glmDF$age<-as.numeric(as.character(glmDF$age))
glmDF$bmi<-as.numeric(as.character(glmDF$bmi))
```

Initialize a matrix to store the results of GLM models<br>

**_Input_**

```r
resGLM<-matrix(nrow=0, ncol=6)
```

Loop through each taxonomic group and fit a GLM with covariates BMI, sex, age, and case/control status<br>

**_Input_**

```r
for (i in colnames(biomeDataS))
{
  model0<- glm(glmDF[,i][[1]] ~ glmDF[,'bmi'][[1]]+glmDF[,'sex'][[1]]+glmDF[,'age'][[1]]+glmDF[,'case_control'][[1]])
  tr<-summary(model0)
  tr<-tr$coefficients
  tr[,'Pr(>|t|)']<-round(as.numeric(tr[,'Pr(>|t|)']),2)
  tr<-cbind(rep(i,nrow(tr)),tr)
  tr<-cbind(c('intercept','BMI','SEX','AGE','control_vs_case'),tr)
  rownames(tr)<-c()
  resGLM<-rbind(resGLM, tr)
}
```

Assign column names to the GLM results and adjust p-values using FDR method<br>

**_Input_**

```r
colnames(resGLM)<-c('factor','tax','estimate','std_error','t_val','pval')
resGLM<-as.data.frame(resGLM)
resGLM<-cbind(resGLM, p.adjust(resGLM$pval, method = 'fdr'))
colnames(resGLM)[length(colnames(resGLM))]<-'pval_adj'
```

Filter results to show only significant associations with case/control status<br>

**_Input_**

```r
resGLM<-resGLM[which(resGLM$factor == 'control_vs_case'),]
resGLMFilt<-resGLM[which(resGLM$pval_adj<0.05),]
resGLMFilt
```

**_Output_**

```
              factor
5    control_vs_case
325  control_vs_case
405  control_vs_case
420  control_vs_case
480  control_vs_case
570  control_vs_case
625  control_vs_case
1055 control_vs_case
                                                                                                     tax
5        Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter
325         Bacteria;Firmicutes;Clostridia;Clostridiales;(Eubacteriaceae/Lachnospiraceae)_2;unclassified
405                              Bacteria;Firmicutes;Clostridia;Clostridiales;Eubacteriaceae;Eubacterium
420        Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;(Lachnoclostridium/unclassified)
480                               Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia
570                                Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Gemmiger
625  Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Candidatus_Stoquefichus
1055        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
                estimate          std_error             t_val pval
5    -0.0472473208760292 0.0136202994412748  -3.4688900254903    0
325     1.49958472903816  0.488411879979953   3.0703281195775    0
405    0.519189420564998  0.181942484855474  2.85359090801369    0
420    0.132428236193973 0.0406417659361785  3.25842721504598    0
480     1.19150475309103  0.332437879526919  3.58414256157155    0
570   0.0341504148648222 0.0118551276702192  2.88064505206554    0
625   -0.287410863591875 0.0917414127348212 -3.13283668764331    0
1055   -2.73511147819197  0.826501709219534 -3.30926294245022    0
     pval_adj
5           0
325         0
405         0
420         0
480         0
570         0
625         0
1055        0
```

Check for overlap between GLM and Wilcoxon test results (significant taxa)<br>

**_Input_**

```r
resGLMFilt$tax[which(resGLMFilt$tax %in% wilcoxResSign$tax)]
```

**_Output_**

```
[1] "Bacteria;Firmicutes;Clostridia;Clostridiales;(Eubacteriaceae/Lachnospiraceae)_2;unclassified" 
[2] "Bacteria;Firmicutes;Clostridia;Clostridiales;Eubacteriaceae;Eubacterium"                      
[3] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;(Lachnoclostridium/unclassified)"
[4] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia"                       
[5] "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia" 
```

## **Step 3: PERMANOVA Analysis**

PERMANOVA to test if there is a significant difference in community composition between cases and controls<br>

Calculate Bray-Curtis distance matrix on biome data<br>

**_Input_**

```r
biomeDist<- DistContinuous(biomeDataS, coef = 'Bray_Curtis')
```

Perform PERMANOVA on the distance matrix using case/control status as the grouping factor<br>

**_Input_**

```r
biomePERM=PERMANOVA(biomeDist, as.factor(metaDataS$case_control),  CoordPrinc = TRUE, PostHoc = 'fdr')
biomePERM$pvalue
```
**_Output_**

```
[1] 0.001998002
```

**_Input_**

```r
summary(biomePERM)
```

**_Output_**

```
 ###### PERMANOVA Analysis #######

Call
PERMANOVA(Distance = biomeDist, group = as.factor(metaDataS$case_control), 
    CoordPrinc = TRUE, PostHoc = "fdr")
________________________________________________

Contrast Matrix
          Case Control
C Case       1       0
C Control    0       1
________________________________________________

PerMANOVA
      Explained Residual df Num df Denom    F-exp     p-value
Total 0.6153953 50.81976      1      283 3.426952 0.001998002
      p-value adj.
Total  0.001998002
________________________________________________

Contrasts
          Explained Residual df Num df Denom    F-exp     p-value
C Case    0.2504767 50.81976      1      283 1.394829 0.001998002
C Control 0.3649186 50.81976      1      283 2.032122 0.001998002
Total     0.6153953 50.81976      1      283 3.426952 0.001998002
          p-value adj.
C Case     0.002000000
C Control  0.002000000
Total      0.001998002
________________________________________________
```

## **Step 4: Balance and Compositional Analysis (`Selbal`)**

Add a small constant (pseudo count) to the biomeDataS to avoid log(0) issues for balance analysis<br>

**_Input_**

```r
biomeDataS_pseudo<-biomeDataS + 0.0001
```

Perform Nearest Balance (NB) analysis to identify microbial taxa associated with case/control status<br>
`abundance` is the transformed microbiome data, `metadata` contains case/control labels<br>

**_Input_**

```r
nb_2 <- nb_lm(abundance = biomeDataS_pseudo,
              metadata = metaDataS,
              pred = "case_control")
```

Retrieve the identified taxa that contribute to the balance in cases and controls<br>

**_Input_**

```r
nb_2$nb$b1$num # Taxa associated with case group
```

**_Output_**

```
 [1] "Bacteria;Firmicutes;Clostridia;Clostridiales;(Eubacteriaceae/Lachnospiraceae)_2;unclassified"          
 [2] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anaerostipes"                             
 [3] "Bacteria;Firmicutes;Clostridia;Clostridiales;(Eubacteriaceae/Lachnospiraceae);unclassified"            
 [4] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;(Lachnoclostridium/unclassified)"         
 [5] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia"                                
 [6] "Bacteria;Firmicutes;(Bacilli/Clostridia);unclassified;unclassified;unclassified"                       
 [7] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;(Blautia/unclassified)"                   
 [8] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospira"                              
 [9] "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;unclassified"                              
[10] "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;Pasteurellaceae;Haemophilus"                
[11] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Gemmiger"                                 
[12] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium"                         
[13] "Bacteria;Firmicutes;Clostridia;Clostridiales;Eubacteriaceae;Eubacterium"                               
[14] "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Serratia"             
[15] "Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Holdemania"                
[16] "Bacteria;Firmicutes;Clostridia;Clostridiales;Peptostreptococcaceae;unclassified"                       
[17] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Micrococcaceae;Rothia"                          
[18] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;unclassified"                             
[19] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Tyzzerella"                               
[20] "Bacteria;Proteobacteria;Alphaproteobacteria;RF32;unclassified;unclassified"                            
[21] "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia"          
[22] "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Comamonas"                   
[23] "Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter"               
[24] "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Alcaligenaceae;Alcaligenes"                 
[25] "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Providencia"          
[26] "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;unclassified"         
[27] "Bacteria;Firmicutes;Clostridia;Clostridiales;(Clostridiaceae/unclassified);unclassified"               
[28] "Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Lactococcus"                              
[29] "Bacteria;Firmicutes;Bacilli;Lactobacillales;Carnobacteriaceae;Granulicatella"                          
[30] "Bacteria;Firmicutes;Clostridia;Clostridiales;Peptococcaceae;Peptococcus"                               
[31] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;(Ruminococcus/unclassified)"              
[32] "Bacteria;Firmicutes;unclassified;unclassified;unclassified;unclassified"                               
[33] "Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella"                          
[34] "Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Xanthomonadaceae;Stenotrophomonas"         
[35] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Blautia"                                  
[36] "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;(Diaphorobacter/Hylemonella)"
[37] "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium"                               
[38] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnoclostridium" 
```

**_Input_**

```r
nb_2$nb$b1$den # Taxa associated with control group
```

**_Output_**

```
 [1] "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia"                                  
 [2] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Porphyromonas"                                             
 [3] "Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella"                                              
 [4] "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales_Family_XI._Incertae_Sedis;Anaerococcus"                             
 [5] "Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;unclassified"                                                 
 [6] "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales_Family_XI._Incertae_Sedis;Peptoniphilus"                            
 [7] "Bacteria;Firmicutes;Clostridia;Clostridiales;(Ruminococcaceae/unclassified);unclassified"                                      
 [8] "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus"                                                    
 [9] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Corynebacteriaceae;Corynebacterium"                                     
[10] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Mobiluncus"                                            
[11] "Bacteria;Proteobacteria;Deltaproteobacteria;Desulfovibrionales;Desulfovibrionaceae;Desulfovibrio"                              
[12] "Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Candidatus_Stoquefichus"                           
[13] "Bacteria;Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium"                                   
[14] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Paraprevotella"                                                
[15] "Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Varibaculum"                                           
[16] "Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter"                               
[17] "Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter"                              
[18] "Bacteria;Proteobacteria;Deltaproteobacteria;Desulfovibrionales;Desulfovibrionaceae;Bilophila"                                  
[19] "Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera"                                                  
[20] "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;unclassified"                                        
[21] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Odoribacter"                                               
[22] "Bacteria;Firmicutes;Clostridia;Clostridiales;[Mogibacteriaceae];unclassified"                                                  
[23] "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales_Family_XI._Incertae_Sedis;(Finegoldia/Unclassified_Tissierellaceae)"
[24] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Butyricimonas"                                             
[25] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Anoxystipes"                                                      
[26] "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus"                                                       
[27] "Bacteria;Firmicutes;Clostridia;Clostridiales;Catabacteriaceae;Catabacter"                                                      
[28] "Bacteria;Synergistetes;Synergistia;Synergistales;Synergistaceae;Cloacibacillus"                                                
[29] "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Proteus"                                      
[30] "Bacteria;Actinobacteria;Actinobacteria;Coriobacteriales;Coriobacteriaceae;Gordonibacter"                                       
[31] "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;[Ruminococcus]"                                                   
[32] "Bacteria;Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Holdemanella"                                      
[33] "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;(Citrobacter/Raoultella)"                     
[34] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;(Faecalibacterium/Gemmiger)"                                      
[35] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides"                                           
[36] "Bacteria;Firmicutes;Negativicutes;Acidaminococcales;Acidaminococcaceae;Acidaminococcus"                                        
[37] "Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Oxalobacteraceae;Oxalobacter"                                       
[38] "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes"                                                      
[39] "Bacteria;Lentisphaerae;Lentisphaeria;Victivallales;Victivallaceae;Victivallis"                                                 
[40] "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Candidatus_Soleaferrea" 
```

Create a heatmap showing the association of microbial taxa with case and control samples<br>

**_Input_**

```r
heatmap_with_split(
  abundance = biomeDataS_pseudo,         # Pseudocount-adjusted biome data
  metadata = metaDataS,                  # Metadata with case/control status
  formula = ~ case_control,              # The grouping factor (case vs control)
  balance = nb_2$nb$b1,                  # Taxa that define the balance
  show_samp_names = FALSE,               # Hide sample names for clarity
  num_name = "taxa_case",                # Label for case-associated taxa
  den_name = "taxa_control",             # Label for control-associated taxa
  others_name = "not associated with Parkinson"  # Other taxa
)
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_03_Parkinsons_disease/heatmap.png" width="100%">

???+ note
    This heatmap is not perfect at all. It is not in "publication ready quality". Yet it serves well for the aim of analysis.

## **Step 5: Machine Learning: Random Forest**

Random Forest classifier for case/control prediction<br>

Split the data into training (80%) and testing (20%) sets<br>

**_Input_**

```r
train_index <- createDataPartition(y = metaDataS$case_control, p = 0.8, list=FALSE)
train_set <- biomeDataS[train_index,]
test_set <- biomeDataS[-train_index,]
```

Train the Random Forest model with 500 trees<br>

**_Input_**

```r
modRF <- randomForest(train_set,as.factor(metaDataS$case_control[train_index]), nTree=500)
```

Evaluate the importance of each variable (taxon) in the Random Forest model<br>

**_Input_**

```r
imp<-as.data.frame(importance(modRF))
imp<-cbind(rownames(imp),imp)
imp<-imp[order(imp$MeanDecreaseGini, decreasing =T),]
```

???+ note
    Run the cell below by yourself to see the dataframe with the importance of each variable (taxon) in the Random Forest model

**_Input_**

```r
View(imp)
```

Plot variable importance (Mean Decrease in Gini index)<br>

**_Input_**

```r
varImpPlot(modRF)
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_03_Parkinsons_disease/varimpplot.png" width="100%">

Predict case/control status for the test set using the trained Random Forest model<br>

**_Input_**

```r
predictionTree <- predict(modRF, newdata=test_set, type='prob')
```

Convert the predicted probabilities into percentages and merge with actual case/control labels<br>

**_Input_**

```r
predictionTreePerc <- as.data.frame(100*predictionTree)
predictionTreePerc<-cbind(rownames(biomeDataS)[-train_index],predictionTreePerc)
predictionTree <- predict(modRF, newdata=test_set)
predictionTreePerc<-cbind(predictionTreePerc,predictionTree)
colnames(predictionTreePerc)[1]<-'sample_name'
```

Merge predictions with the actual case/control status from metadata<br>

**_Input_**

```r
predictionTreePerc<-inner_join(predictionTreePerc, metaDataS[,c('sample_name','case_control')])
```

Confusion matrix to assess Random Forest performance<br>

**_Input_**

```r
confmatRF<-confusionMatrix(data=as.factor(predictionTreePerc$predictionTree), reference=as.factor(predictionTreePerc$case_control))
confmatRF
```

**_Output_**

```
Confusion Matrix and Statistics

          Reference
Prediction Case Control
   Case      25      12
   Control    8      11
                                          
               Accuracy : 0.6429          
                 95% CI : (0.5036, 0.7664)
    No Information Rate : 0.5893          
    P-Value [Acc > NIR] : 0.2501          
                                          
                  Kappa : 0.2422          
                                          
 Mcnemar's Test P-Value : 0.5023          
                                          
            Sensitivity : 0.7576          
            Specificity : 0.4783          
         Pos Pred Value : 0.6757          
         Neg Pred Value : 0.5789          
             Prevalence : 0.5893          
         Detection Rate : 0.4464          
   Detection Prevalence : 0.6607          
      Balanced Accuracy : 0.6179          
                                          
       'Positive' Class : Case     
```

### **Make ROC (Receiver Operating Characteristic) curve and calculate AUC (Area Under Curve)**

Obtain predicted probabilities for the test set for ROC curve creation<br>

**_Input_**

```r
prediction_for_roc_curve <- predict(modRF, newdata=test_set, type='prob')
classes <- levels(as.factor(predictionTreePerc$case_control))
```

Use the second column of predicted probabilities (i.e., for the "case" class)<br>

**_Input_**

```r
pred <- prediction(prediction_for_roc_curve[,2],predictionTreePerc$case_control)
```

Create a performance object for the true positive rate (TPR) and false positive rate (FPR)<br>

**_Input_**

```r
perf <- performance(pred, "tpr", "fpr")
```

Calculate the AUC for the ROC curve<br>

**_Input_**

```r
aucRF <- performance(pred, "auc")
aucRF <- aucRF@y.values[[1]]
aucRF = round(aucRF,2)
```

Plot the ROC curve with AUC displayed<br>

**_Input_**

```r
plot(perf,main="ROC Curve", sub = paste("AUC = ",aucRF) )
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_03_Parkinsons_disease/roccurve.png" width="100%">