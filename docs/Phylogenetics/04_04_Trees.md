# **Preparing the alignment and building trees**

## **Instruction**
![CLT](https://img.shields.io/badge/Language-Command Line Tools & Python-steelblue)<br>
![IDE](https://img.shields.io/badge/Recommended IDE-Jupyter Notebook-steelblue)

???+ question
    Where to get the data?

Download it from GitHub repository:

```bash
wget https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/data/04_Phylogenetics/04_04_Trees.zip
```

```bash
unzip 04_04_Trees.zip && rm -rf 04_04_Trees.zip
```

These are the data we will be working with:

```
Archive:  04_04_Trees.zip
   creating: data/
  inflating: data/SUP35_10seqs.fa
```

>For this work, we will use the alignment of [SUP35 gene](https://www.yeastgenome.org/locus/S000002579) obtained by `prank` considering codons.<br>

**_Input_**

```bash
prank -codon -d=data/SUP35_10seqs.fa -o=data/SUP35_aln_prank.best.fas -F
```

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_04_Trees/tree_building_overview.png" align='center', width="100%">
</div>

_Tree building overview_

----------------------------------------------

## **Step 1: Cut bad areas out of the alignment using `trimAl`**

**_Input_**

```bash
trimal -in data/SUP35_aln_prank.best.fas -out data/SUP35_aln_prank.trim.fas -automated1
```

----------------------------------------------

## **Step 2: fit an evolution model in `ModelTest` (`ModelTest-NG`)**

**_Input_**

```bash
modeltest-ng -i data/SUP35_aln_prank.trim.fas -o data/modeltest/SUP35_trim_modeltest
```

**_Output_**

```
                     Model         Score        Weight
   BIC             TIM3+G4    18180.5614        0.3950
   AIC           TIM3+I+G4    18041.1550        0.5377
  AICc           TIM3+I+G4    18041.1550        0.5377
```

In total we see that the `TIM3+G4` is recognised as the best model!

----------------------------------------------

## **Step 3: Build an ML-tree in RAxML-NG using the selected model**

First, let's check that our tree is being built at all

**_Input_**

```bash
raxml-ng --check --msa data/SUP35_aln_prank.trim.fas  --model TIM3+G4 --prefix data/raxml/SUP35_raxml_test
```

**_Output_**

```
Alignment can be successfully read by RAxML-NG.
```

Great! Let's go!

**_Input_**

```bash
raxml-ng --msa data/SUP35_aln_prank.trim.fas --model TIM3+G4 --prefix data/raxml/SUP35_raxml --threads 2 --seed 222  --outgroup SUP35_Kla_AB039749
```

----------------------------------------------

## **Step 4: Draw the resulting tree (the best ML tree)**

For drawing trees we will use a simple prewritten R script. Let's take a look at it to understand how it works:

``` R title="draw_tree.R"
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(ggtree)
tr <- read.tree(args[1]) ##SUP35_raxml.raxml.bestTree
ggtree(tr) + geom_tiplab() + xlim(0,2) + 
  geom_treescale()

ggsave(args[2])
```

**_Input_**

```bash
Rscript scripts/draw_tree.R data/raxml/SUP35_raxml.raxml.bestTree imgs/SUP35_raxml.png
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_04_Trees/SUP35_raxml.png" align='center', width="100%">
</div>

----------------------------------------------

## **Step 5: Model selection in ModelFinder (IQ-TREE)**

**_Input_**

```bash
iqtree2 -m MFP -s data/SUP35_aln_prank.trim.fas --prefix data/modelfinder/SUP35_MF2
head -42 data/modelfinder/SUP35_MF2.iqtree | tail -6
```

**_Output_**

```
Best-fit model according to BIC: TIM3+F+G4

List of models sorted by BIC scores: 

Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC
TIM3+F+G4         -8993.686   18035.372 +   0.0517   18035.972 +   0.0549   18170.092 +    0.737
```

We see that the model `TIM3+F+G4` is recognised as the best!

???+ question
    Do the models selected by ModelTest and ModelFinder differ, and how much?

???+ success
    In general, we got the same thing. Only ModelFinder also threw in information about the empirical frequencies of the letters themselves in the alignment.

    |    |ModelTest|ModelFinder|
    |----|---------|-----------|
    |Model|TIM3+G4|TIM3+F+G4|
    |BIC|18180.5614|18170.092|

----------------------------------------------

## **Step 6: Build an ML tree in IQ-TREE using the selected model**

**_Input_**

```bash
iqtree2 -m TIM3+F+G4 -s data/SUP35_aln_prank.trim.fas --prefix data/iqtree/SUP35_iqtree
```

----------------------------------------------

## **Step 7: Draw the resulting tree (the best ML tree)**

**_Input_**

```bash
Rscript draw_tree.R data/iqtree/SUP35_iqtree.treefile imgs/SUP35_iqtree.png
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_04_Trees/SUP35_iqtree.png" align='center', width="100%">
</div>

The trees obtained with RAxML and IQTREE are fundamentally similar.<br>
They have different views on how well the outer groups diverge, and they have different topologies.

----------------------------------------------

## **Step 8: Comparison of likelihood (log likelihood) of trees obtained with different models and before and after filtering**

???+ question
    What conclusion can be drawn from this?

**_Input_**

```bash
iqtree2 -s data/SUP35_aln_prank.best.fas -pre data/iqtree_unfilt/SUP35_iqtree_unfilt
grep "Log-likelihood" data/iqtree_unfilt/SUP35_iqtree_unfilt.iqtree
```

**_Output_**

```
Log-likelihood of the tree: -9696.9044 (s.e. 160.3706)
```

**_Input_**

```bash
grep "Log-likelihood" data/iqtree/SUP35_iqtree.iqtree
```

**_Output_**

```
Log-likelihood of the tree: -8993.1633 (s.e. 149.1347)
```

**_Input_**

```bash
iqtree2 -s data/SUP35_aln_prank.best.fas -m JC -pre data/iqtree_JC/SUP35_iqtree_JC
grep "Log-likelihood" data/iqtree_JC/SUP35_iqtree_JC.iqtree
```

**_Output_**

```
Log-likelihood of the tree: -10482.7253 (s.e. 176.1729)
```

|              |Unfilt|TIM3+F+G4|JC|
|--------------|------|---------|--|
|Log-likelihood|-9696.9044|-8993.1633|-10482.7253|

- Before filtering Log-likelihood = -9696.9044
- When using the `TIM3+F+G4` model = -8993.1633
- When using `JC` model = -10482.7253<br>

???+ success
    That said, the topology is the same everywhere! Anyway, if we only need to look at the topology, we can run iqtree with either model...

----------------------------------------------

## **Step 9: Generation of 100 replicas of a regular bootstrap**

**_Input_**

```bash
time iqtree2 -s data/SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre data/iqtree_bootstrap/SUP35_TIM3_b -b 100
```

----------------------------------------------

## **Step 10: Generation of 1000 ultrafast bootstrap replicas**

**_Input_**

```bash
time iqtree2 -s data/SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre data/iqtree_ultrafast_bootstrap/SUP35_TIM3_ufb -bb 1000
```

???+ question
    What is the difference between normal and ultrafast bootstrap runtimes and the values obtained?

**_Input_**

```python
bootstrap_100 = 3 * 60 + 17
ultrafast_bootstrap_1000 = 3.258

print(f'Generation of ultrafast bootstrap is faster: {bootstrap_100 / ultrafast_bootstrap_1000} times')
```

**_Output_**

```
Generation of ultrafast bootstrap is faster: 60.46654389195826 times
```

???+ success
    Generating 100 replicas of regular bootstrap took 3:17.00, while generating 1000 replicas of ultrafast bootstrap took 3.258. That's a huge difference!

----------------------------------------------

## **Step 11: Running the previous command, but with generation: 1000 ultrafast bootstrap + 1000 alrt + abayes**

**_Input_**

```bash
iqtree2 -s data/SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre data/iqtree_ufb_alrt_abayes/SUP35_TIM3_B_alrt_abayes -bb 1000 -alrt 1000 -abayes
```

----------------------------------------------

## **Step 12: Drawing the resulting tree with three supports**

For this step we will use a bit different R script. Let's take a look at it:

``` R title="draw_tree_max.R"
args <- commandArgs(trailingOnly=TRUE)

library(ggtree)
tree_alrt_abayes_ufb <- read.tree(args[1])
ggtree(tree_alrt_abayes_ufb) + 
  geom_tiplab() + geom_nodelab() +
  geom_treescale() + xlim(0, 0.9)
# funny labels
label <- tree_alrt_abayes_ufb$node.label
alrt <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 1)) #sun
abayes <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 2)) #yin yang
ufb <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 3)) #star
large_alrt <- ifelse(alrt > 70, intToUtf8(9728), "")
large_abayes <- ifelse(abayes > 0.7, intToUtf8(9775), "")
large_ufb <- ifelse(ufb > 95, intToUtf8(9733), "")
newlabel <- paste0(large_alrt, large_abayes, large_ufb)
tree_alrt_abayes_ufb$node.label <- newlabel
ggtree(tree_alrt_abayes_ufb) + 
  geom_tiplab() + geom_nodelab(nudge_x = -.01, nudge_y = .1) +
  geom_treescale() + xlim(0, 0.9)

ggsave(args[2])
```

**_Input_**

```bash
Rscript scripts/draw_tree_max.R data/iqtree_ufb_alrt_abayes/SUP35_TIM3_B_alrt_abayes.treefile imgs/SUP35_TIM3_B_alrt_abayes.png
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_04_Trees/SUP35_TIM3_B_alrt_abayes.png" align='center', width="100%">
</div>

All values ​​- `alrt`, `abayes` and `ufb` are positively correlated - that is, if the indicators are high, then they are all high, as a rule. But this is not proportional.<br>
In the tree above, the indicator values ​​are replaced with the symbols of the sun, yin-yang and asterisk for educational purposes. In real life they don’t use this, but giving bare numbers through a slash is also not comme il faut...<br>