# **Practice 1: Crohn's disease**
>This chapter contains a training data for 16S amplicon analysis and an example of analysis.

## **Introduction**

In this study, the gut microbiome composition in individuals with Crohn’s Disease (CD) and Healthy Controls (HC) was investigated. The goal was to identify taxonomic shifts, potential functional implications, and associations with disease. High-throughput sequencing data and various bioinformatics tools were employed to analyze microbial diversity and abundance.<br>

Data obtained from the article [«_A microbial signature for Crohn's disease_»](https://pubmed.ncbi.nlm.nih.gov/28179361/) was used. These data include information on microbial samples from healthy individuals and patients with Crohn's disease.<br>

## **Instruction**

You can run commands below in your `RStudio` in `R script`.<br>
Or if you want to write a beautiful & convenient to read laboratory journal you can use `R Markdown`.<br>

----------------------------------------------

## **Step 1: Loading libraries and data**

In this 16S amplicon analysis pipeline we will use `MicrobeR`, `balance`, `NearestBalance` & `selbal` packages. Their installation is a bit difficult. Please follow the code below:<br>

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("philr")
BiocManager::install("DECIPHER")

devtools::install_github("jbisanz/MicrobeR")
devtools::install_github("tpq/balance")
devtools::install_bitbucket("knomics/nearestbalance")
devtools::install_github(repo = "malucalle/selbal")
```

First, load all the libraries needed for the data analysis.<br>

**_Input_**

```r
library(data.table)
library(openxlsx)
library(MicrobeR)
library(ggplot2)
library(zCompositions)
library(NearestBalance)
library(GUniFrac)
library(vegan)
library(ape)
library(selbal)
```
 
Then set the working directory.<br>
 
**_Input_**

```r
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)
```

Download the data to work with.<br>

**_Input_**

```r
url <- "https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/data/05_16S_amplicon_analysis/05_02_Crohns_disease.zip"

zipF<- "05_02_Crohns_disease.zip"

download.file(url, zipF)

outDir<-"."

unzip(zipF,exdir=outDir)

if (file.exists(zipF)) {
  file.remove(zipF)
}
```

Load metadata and sort it by participant number.<br>

**_Input_**

```r
metadata = fread("data/metadata.csv")
metadata[, .N, by = diagnosis_full]
```

**_Output_**

```
diagnosis_full  N
<chr> <int>
CD  34
HC  34
2 rows
```

Load the `Counts` table.<br>

**_Input_**

```r
counts <- read.csv("data/counts.csv",row.names = 1)
counts <- counts[metadata$sample, colSums(counts) >0]
```

----------------------------------------------

## **Step 2: Check the data**

How many samples and microbial samples?<br>

**_Input_**

```r
dim(counts)
```

**_Output_**

```
[1]  68 210
```

What's the coverage?<br>

**_Input_**

```r
range(rowSums(counts))
```

**_Output_**

```
[1]  19414 121234
```

Composition of samples.<br>

**_Input_**

```r
raw_abund <- counts/rowSums(counts)
```

**_Input_**

```r
metadata$SampleID <- metadata$sample
```

**_Input_**

```r
Microbiome.Barplot(t(raw_abund), metadata, CATEGORY = "diagnosis_full")
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/microbiome_barplot.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/microbiome_barplot.jpg", width = 11, height = 3.5, dpi = 300)
```

Is there enough coverage?<br>

**_Input_**

```r
metadata$coverage <- rowSums(counts)
```

**_Input_**

```r
ggplot(metadata) + 
  geom_histogram(aes(coverage)) + 
  theme_minimal() + 
  xlab("N samples")
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/coverage_quality.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/coverage_quality.jpg", width = 11, height = 3.5, dpi = 300)
```

Sequencing quality is sufficient.<br>

----------------------------------------------

## **Step 3: Filtration from rare and under-represented taxa**

Keeping the microbes that occur in >30% of samples.<br>

**_Input_**

```r
filt_counts <- counts[, colSums(counts>0)>0.3*nrow(counts)]
metadata[, filt_coverage := rowSums(filt_counts)]
metadata[, proportion_of_prevalent_taxa := 100*filt_coverage/coverage]
```

Coverage of samples after filtration.<br>

**_Input_**

```r
ggplot(metadata)+
  geom_histogram(aes(filt_coverage)) + 
  theme_minimal()+
  xlab("Post-filtration coverage") + 
  ylab("N samples")
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/post-filtration_coverage.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/post-filtration_coverage.jpg", width = 11, height = 3.5, dpi = 300)
```

What proportion of the microbes remained in the analysis.<br>

**_Input_**

```r
ggplot(metadata)+
  geom_histogram(aes(proportion_of_prevalent_taxa)) + 
  theme_minimal()+ 
  xlab("Proportion of microbes remaining in the assay") + 
  ylab("N samples")
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/remaining_proportion.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/remaining_proportion.jpg", width = 11, height = 3.5, dpi = 300)
```

How many samples and microbial samples?<br>

**_Input_**

```r
dim(filt_counts)
```

**_Output_**

```
[1] 68 89
```

What's the coverage?<br>

**_Input_**

```r
range(rowSums(filt_counts))
```

**_Output_**

```
[1]  12981 120882
```

**_Input_**

```r
matrix_data <- matrix(c("Before", 210, 19412,
                         "After", 89, 12981), 
                      nrow = 2, byrow = TRUE)

data <- as.data.frame(matrix_data)

colnames(data) <- c("", "N microbes", "Minimum coverage")

print(data)
```

**_Output_**

```
  N microbes  Minimum coverage
<chr> <chr> <chr>
Before	210	19412
After	89	12981
```

Calculating relative abundances.<br>

**_Input_**

```r
abundance <- cmultRepl(filt_counts)
```

**_Output_**

```
No. adjusted imputations:  694 
```

**_Input_**

```r
heatmap_with_split(abundance,
                   metadata, ~ diagnosis_full,
                   show_samp_names = F) + 
  theme(axis.text.y = element_text(size =5))
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/relative_abundance.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/relative_abundance.jpg", width = 8, height = 10, dpi = 300)
```

----------------------------------------------

## **Step 4: Counting alpha diversity**

**_Input_**

```r
alpha_div <- rowMeans(sapply(1:5, function(i){
  counts_rar_i = Rarefy(counts, 19000)$otu.tab.rff
  alpha_div_i = vegan::diversity(counts_rar_i)
}))
metadata$Shannon.index <- alpha_div[metadata$sample]
```

**_Input_**

```r
ggplot(metadata) + 
  geom_boxplot(aes(diagnosis_full, Shannon.index, fill = diagnosis_full)) + 
  theme_minimal() +
  theme(legend.position = 'none') + 
  xlab("")
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/alpha_diversity.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/alpha_diversity.jpg", width = 8, height = 8, dpi = 300)
```

Is it different?<br>

**_Input_**

```r
wilcox.test(Shannon.index ~ diagnosis_full, metadata)$p.value
```

**_Output_**

```
[1] 1.314773e-09
```
The _p_-value of `1.314773e-09` indicates a significant difference in alpha diversity between the two groups. CD samples exhibit lower alpha diversity compared to HC samples.<br>

----------------------------------------------

## **Step 5: Aitchison's beta diversity: is there a difference in proportions?**

**_Input_**

```r
clr <- log(abundance) - rowMeans(log(abundance))
beta_div <- dist(clr)
```

**_Input_**

```r
pcoa_res <- pcoa(beta_div)$vectors
var <- apply(pcoa_res, 2, var)
var_rel <- round(var*100/sum(var), 1)
```

**_Input_**

```r
ggplot(cbind(metadata,pcoa_res)) + 
  geom_point(aes(Axis.1, Axis.2, col=diagnosis_full)) +
  coord_fixed() + 
  theme_minimal() + 
  labs(col="") + 
  xlab(paste0("Axis.1 (",var_rel[1], "%)")) + 
  ylab(paste0("Axis.2 (",var_rel[2], "%)"))
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/beta_diversity.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/beta_diversity.jpg", width = 8, height = 8, dpi = 300)
```

Is the difference statistically significant?<br>

**_Input_**

```r
adonis2(beta_div ~ metadata$diagnosis_full)
```

**_Output_**

```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = beta_div ~ metadata$diagnosis_full)
                        Df SumOfSqs      R2      F Pr(>F)    
metadata$diagnosis_full  1     5400 0.16562 13.101  0.001 ***
Residual                66    27205 0.83438                  
Total                   67    32605 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

The statistical analysis using the `adonis2` test shows that the difference in beta diversity between the two groups is statistically significant (p-value = 0.001).<br>

----------------------------------------------

## **Step 6: What exactly is the difference?**

**_Input_**

```r
nb <- nb_lm(abundance = abundance,
            metadata = metadata,
            pred = "diagnosis_full")
```

**_Input_**

```r
heatmap_with_split(abundance = abundance[, unlist(nb$nb$b1)],
                   metadata = metadata,
                   formula = ~ diagnosis_full,
                   num_name = "health-related",
                   den_name = "disease-related",
                   show_samp_names = F,
                   balance = nb$nb$b1)
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/heatmap_w_split.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/heatmap_w_split.jpg", width = 10, height = 10, dpi = 300)
```

Illustrating the difference between the average microbiota of healthy and sick people.<br>

**_Input_**

```r
psi <- make_psi_from_sbp(nb$coord$sbp)
mean_diff_clr <- drop(nb$lm_res$coefficients[2,] %*% psi)
bal_unit <- balance_to_clr(nb$nb$b1, colnames(abundance))
bal_diff_clr <- drop(bal_unit %*% mean_diff_clr) * bal_unit
tab <- data.table(taxon = names(mean_diff_clr),
                  clr_diff = mean_diff_clr,
                  bal = bal_diff_clr)
setorderv(tab, "clr_diff")
tab$taxon <- factor(tab$taxon, levels = tab$taxon)
```

Mean difference between sick and healthy individuals.<br>

**_Input_**

```r
ggplot(tab) + 
  geom_col(aes(clr_diff, taxon), fill = "darkblue") + 
  theme_minimal() + xlab("CLR(v)") + ylab("") + 
  theme(axis.text.y = element_text(size =5))
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/mean_difference.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/mean_difference.jpg", width = 10, height = 10, dpi = 300)
```

Approximate, simplified difference.<br>

**_Input_**

```r
ggplot(tab) + 
  geom_col(aes(bal, taxon), fill = "darkblue") + 
  theme_minimal() + xlab("CLR(b)") + ylab("") +
  theme(axis.text.y = element_text(size =5))
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/approximate_difference.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/approximate_difference.jpg", width = 10, height = 10, dpi = 300)
```

Balance value in each sample.<br>

**_Input_**

```r
metadata$balance <- balance.fromSBP(abundance, nb$nb$sbp)
```

**_Input_**

```r
ggplot(metadata) + 
  geom_boxplot(aes(diagnosis_full, balance, fill = diagnosis_full)) + 
  theme_minimal() +
  theme(legend.position = 'none') + 
  xlab("") 
```

**_Output_**

<img src="/images/ngs-handbook/05_16S_Amplicon_Analysis/05_02_Crohns_disease/balance_value.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/balance_value.jpg", width = 8, height = 8, dpi = 300)
```

Has it changed?<br>

**_Input_**

```r
wilcox.test(balance ~ diagnosis_full, metadata)$p.value
```

**_Output_**

```
[1] 2.621864e-17
```

The p-value of `2.621864e-17` indicates a significant difference in balance between the two groups.<br>
