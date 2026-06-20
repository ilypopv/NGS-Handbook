---
hide:
  - toc
icon: lucide/tree-deciduous
---
# **Phylogenetics**

>This chapter contains a manual on phylogenetic analysis 

## **Instruction**
![language](https://img.shields.io/badge/Language-Command Line Tools, Python & R-steelblue)<br>
![IDE](https://img.shields.io/badge/IDE-Jupyter Notebook & RStudio-steelblue)

To recreate any of the part of this manual please install:

```bash
wget https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/envs/qc.yaml
```

```bash
conda env create -f phylo.yaml
```

And of cource do not forget to activate the envinronment!

```bash
conda activate phylo
```

???+ warning
    Disclaimer: when using `ete3` on `WSL` it will crash with error if you try to render any kind of picture (unless you have tunned `X11` forwarding). It is recommended to use `Ubuntu` or `macOS` for steps with `ete3`. Anyway there is always an alternative!

----------------------------------------------

## **08 Tanglegram**

In the [Tanglegram chapter](04_08_Tanglegram) there is guide how to make co-phylogeny plots using `R` by several ways: (1) manually creating it with `ggtree`; (2) easily creating small tanglegram with `cophyloplot()` function in `ape` and (3) creating large tanglegram with `dendextend`.

----------------------------------------------

## **07 Visualization Pro**

In the [Visualization Pro chapter](04_07_Visualization_Pro) there is guide how to make a 1st class publication ready quality tree visualization using `ggtree` in `R`. By the end of the guide the circular ML-tree annotated with `tip-points`, `tip-labels` and `heatmaps` is created. The ideas from the first lesson continue here.

----------------------------------------------

## **06 Phylogenomics**

In the [Phylogenomics chapter](04_06_Phylogenomics) there is guide how to reconstruct phylogeny using more than just one gene with two examples of pipelines. The 1st is about using mitochondrial genes and manual work with data. The 2nd is about using proteomes and extracting single-copy orthologs with `Proteinortho`. The main concern in phylogenomics studies is storing the data and naming the sequences properly. This manual pays a lot attention to it and provides insights on using `BBMap` to achieve this goal.

----------------------------------------------

## **05 Rooting and comparing trees & Dating**

In the [Root Date chapter](04_05_Root_Date) there is guide how to root trees using different approaches (by a known external clade / `midpoint rooting` / `non-reversible` model), how to compare these approaches and how to perform root-supporte tree visualisation (`rootstrap`). Also there is a mini-project on dating the common ancestor of smoky leopards with "full-circle" study: downloading the data from NCBI using `efetch`, aligning the sequences with `mafft`, trimming with `trimal` and tree construnction with `iqtree2` followed by visualization and analysis in GUI apps: `Beauti`, `Tracer`, `TreeAnnotator` and `FigTree`.

----------------------------------------------

## **04 Preparing the alignment and building trees**

In the [Trees chapter](04_04_Trees) there is guide how to prepare the alignment by cutting bad areas out of with `trimAl`, how to select the model with `ModelTest-NG` and `ModelFinder`, how to perform a tree search with `RaxML-NG` and `IQ-TREE` and how to do topology assesment with `bootstrap`.

----------------------------------------------

## **03 Multiple sequence alignment**

In the [MSA chapter](04_03_MSA) there is guide how to perform multiple sequence alignment using all possible tools (`clustalw`, `clustalo` `muscle`, `mafft`, `kalign`, `tcoffee`, `prank`).

----------------------------------------------

## **02 Working with NCBI**

In the [Mining data chapter](04_02_Mining_Data) there is guide how to work with NCBI using all possible ways - `R` (`reutils`), `Python` (`Bio::Entrez`) and `bash` (`esearch`).

----------------------------------------------

## **01 Drawing Trees**

In the [Intro Trees chapter](04_01_Intro_Trees) there is guide how to draw phylogenetic trees using both `R` (`ape`, `ggtree`) and `Python` (`Bio::Phylo`, `ete3`).
