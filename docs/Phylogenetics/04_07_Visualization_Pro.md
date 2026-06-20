# **Visualization Pro**

>For this work, we will use sequences and metadata from the viral hemorrhagic septicemia virus (VHSV). This virus is a fish novirhabdovirus (negative stranded RNA virus) with an unusually broad host spectra: it has been isolated from more than 80 fish species in locations around the Northern hemisphere.<br>
>
>This guide was inspired by: https://github.com/acarafat/tutorials. Thanks for the data! Yet there are some new adjustments!
>
>For more guides on pro visualization please visit [bioconnector's `ggtree` workshop](https://bioconnector.github.io/workshops/r-ggtree.html)

## **Instruction**

This guide is about making the 1st class publication ready quality tree visualization.<br>
It is divided into 2 parts:<br>

1. **Making the tree**<br>
  ![CLT](https://img.shields.io/badge/Language-Command Line Tools & Python-steelblue)<br>
  ![IDE](https://img.shields.io/badge/Recommended IDE-Jupyter Notebook-steelblue)
2. **Visualizing the tree**<br>
  ![CLT](https://img.shields.io/badge/Language-R-steelblue)<br>
  ![IDE](https://img.shields.io/badge/Recommended IDE-RStudio-steelblue)

???+ question
    In this chapter we will use a lot of prewritten scripts and even some images for visualization. Where to get them?

Download it from GitHub repository:

```bash
wget https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/data/04_Phylogenetics/04_07_Visualization_Pro.zip
```

```bash
unzip 04_07_Visualization_Pro.zip && rm -rf 04_07_Visualization_Pro.zip
```

These are the scripts and images we will be working with:

```
Archive:  04_07_Visualization_Pro.zip
   creating: data/
  inflating: data/metadata.csv       
  inflating: data/vhsv.fasta
```

----------------------------------------------

## **Part 1: Making the tree**

### **Step 1: Multiple Sequences Alignment**

Run multiple sequences alignment with `mafft`.<br>

**_Input_**

```bash
mafft data/vhsv.fasta > data/vhsv_mafft.fa
```

----------------------------------------------

### **Step 2: Selecting a model in ModelFinder (IQ-TREE)**

Create a directory for `ModelFinder` output.<br>

**_Input_**

```bash
mkdir data/modelfinder
```

Run `ModelFinder`.<br>

**_Input_**

```bash
iqtree2 -m MFP -s data/vhsv_mafft.fa --prefix data/modelfinder/vhsv_MF2
```

Examine the best model of evolution that is the most suitable for our alignment.<br>

**_Input_**

```bash
head -42 data/modelfinder/vhsv_MF2.iqtree | tail -6
```

**_Output_**

```
Best-fit model according to BIC: TVMe+I+R2

List of models sorted by BIC scores: 

Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC
TVMe+I+R2         -6984.748   14221.497 -  0.00151   14244.406 -   0.0031   14892.963 +    0.719
```

In total we see that the model `TVMe+I+R2` is recognised as the best!<br>

----------------------------------------------

### **Step 3: Build an ML-tree in IQ-TREE using the selected model**

Create a directory for `IQ-TREE` output.<br>

**_Input_**

```bash
mkdir data/iqtree
```

Run `IQ-TREE` with 1000 replicates of `bootstrap` and `alrt`.<br>

**_Input_**

```bash
iqtree2 -s data/vhsv_mafft.fa -m TVMe+I+R2 -pre data/iqtree/vhsv -bb 1000 -alrt 1000
```

Examine the `IQ-TREE` output folder.<br>

**_Input_**

```bash
ls data/iqtree
```

**_Output_**

```
vhsv.bionj   vhsv.contree  vhsv.log	vhsv.splits.nex
vhsv.ckp.gz  vhsv.iqtree   vhsv.mldist	vhsv.treefile
```

Examine `vhsv.treefile` file.<br>

**_Input_**

```bash
cat data/iqtree/vhsv.treefile
```

**_Output_**

```
(AU-8-95:0.0110908144,(CH-FI262BFH:0.0131626927,(DK-200098:0.0013245216,DK-9995144:0.0000021124)99.9/100:0.0095461686)9/64:0.0007239176,((((((((DK-1p40:0.0019819713,DK-1p86:0.0019797009)76.6/97:0.0006555511,((DK-1p8:0.0000009918,(((((DK-4p37:0.0000009918,SE-SVA14:0.0006553525)86.8/99:0.0006552139,DK-5e59:0.0013137944)0/59:0.0000009918,SE-SVA-1033:0.0000009918)78.3/95:0.0006575436,((DK-6p403:0.0006555940,SE-SVA31:0.0006553577)0/84:0.0000009918,UK-MLA98-6HE1:0.0000009918)85.8/98:0.0019746439)74.9/96:0.0006544639,KRRV9601:0.0006553263)0/5:0.0000009918)0/35:0.0000009918,UK-9643:0.0033018138)77.4/95:0.0006559696)99.4/100:0.0054050544,DK-M.rhabdo:0.0026557033)97.1/99:0.0043768013,((((((DK-1p53:0.0006554948,DK-1p55:0.0000009918)100/100:0.0816239127,(US-Makah:0.0285256070,US-Goby1-5:0.0140617136)100/100:0.1210196539)45.4/78:0.0178122396,(((DK-4p101:0.0098266846,((DK-4p168:0.0013124829,(UK-H17-2-95:0.0026484120,(UK-H17-5-93:0.0006595661,UK-MLA98-6PT11:0.0013182714)76.2/99:0.0006379007)76.8/99:0.0006676816)79.4/100:0.0010717553,IR-F13.02.97:0.0055710520)87.9/96:0.0020014929)75.6/98:0.0009426132,FR-L59X:0.0128796948)94.9/100:0.0075738902,UK-860-94:0.0158845820)100/100:0.0465027783)99.6/100:0.0335713781,GE-1.2:0.0179796623)90.6/90:0.0051656186,(((DK-2835:0.0040058973,(DK-5123:0.0028538659,DK-5131:0.0044317727)78/92:0.0031564427)99.9/100:0.0132939283,DK-Hededam:0.0099438696)73.4/75:0.0008365002,DK-F1:0.0114614651)0/60:0.0000024805)79.3/82:0.0007471327,((FI-ka422:0.0032972674,FI-ka66:0.0000009918)92.8/100:0.0026339174,NO-A16368G:0.0019943447)92.3/97:0.0020403429)61.6/78:0.0002733446)98.6/96:0.0054285225,(FR-1458:0.0099883404,FR-2375:0.0066740966)90.8/95:0.0027197676)93.6/99:0.0031601931,(((((((((((DK-200027-3:0.0020089775,DK-200079-1:0.0020149815)76.7/89:0.0006222640,(DK-9795568:0.0019716416,DK-9995007:0.0013117492)0/86:0.0000009918)83/84:0.0006608276,DK-9895093:0.0013141923)0/73:0.0000020835,DK-7380:0.0026333937)52.9/85:0.0012862237,DK-9595168:0.0020301147)95.9/99:0.0026894771,DK-6045:0.0000009918)97.5/100:0.0026369103,DK-7974:0.0039639852)0/85:0.0000009918,DK-5151:0.0000009918)75.6/97:0.0006726697,DK-6137:0.0006450449)99.9/100:0.0115741245,((DK-3592B:0.0019729337,(DK-3946:0.0006552354,Fil3:0.0006618033)0/65:0.0000009918)0/86:0.0000742921,DK-3971:0.0025543911)99.5/100:0.0068250583)0/86:0.0002859106,(DK-5741:0.0000009918,(DK-9695377:0.0013153923,(DK-9895024:0.0000009918,DK-9895174:0.0013114345)85.7/98:0.0013140555)90.9/99:0.0013147689)99.7/100:0.0081963781)95.9/77:0.0050968140)79.8/66:0.0033929129,FR-0284:0.0121797246)78.3/89:0.0028127442,FR-0771:0.0001556512)99.4/100:0.0075113449);
```

This is the tree in `Newick` format that we will use for visualization!<br>

----------------------------------------------

### **Step 4: Examine metadata**

Import `pandas`.<br>

**_Input_**

```python
import pandas as pd
```

Read the metadata file.<br>

**_Input_**

```python
pd.read_csv('data/metadata.csv')
```

**_Output_**

| |Strain|Host|Water|Country|ACCNo|Year|
|-|-|-|-|-|-|-|
|0|AU-8-95|Rainbow trout|Fresh water|AU|AY546570.1|1995|
|1|CH-FI262BFH|Rainbow trout|Fresh water|CH|AY546571.1|1999|
|2|DK-1p40|Rockling|Sea water|DK|AY546575.1|1996|
|3|DK-1p53|Atlantic Herring|Sea water|DK|AY546577.1|1996|
|4|DK-1p55|Sprat|Sea water|DK|AY546578.1|1996|
|...|...|...|...|...|...|...|
|56|UK-H17-5-93|Cod|Sea water|UK|AY546630.1|1993|
|57|UK-MLA98-6HE1|Herring|Sea water|UK|AY546631.1|1998|
|58|UK-MLA98-6PT11|Norway prout|Sea water|UK|AY546632.1|1998|
|59|US-Makah|Coho salmon|Fresh water|USA|U28747.1|1988|
|60|US-Goby1-5|Round goby|Unknown|USA|AB672615.1|2006|

61 rows × 6 columns

----------------------------------------------

## **Part 2: Visualizing the tree**

### **Step 1: Preparation**

1st step - install or call libraries.<br>

**_Input_**

```r
if (!require("pacman")) install.packages("pacman")

pacman::p_load(ggplot2, ggtree, phangorn, treeio, ggnewscale, viridis)
```

2nd step - set the working directory.<br>

**_Input_**

```r
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)
```

----------------------------------------------

### **Step 2: Initial visualization**

Create a directory for images.<br>

**_Input_**

```r
dir.create("imgs", showWarnings = TRUE, recursive = FALSE, mode = "0777")
```

As we remember from the previous laboratory journal our tree file is in `Newick` format and it was created with `IQ-TREE`. Also, there are bootstrap values in our tree files. That is why it is better to read the tree with `read.iqtree` function instead of default `read.tree` function.<br>

So let us read the tree, midpoint root it and make the 1st visualization.<br>

**_Input_**

```r
#Read the tree
vshv.tree <- read.iqtree('data/iqtree/vhsv.treefile')

#Midpoint root the tree
vshv.tree@phylo <- midpoint(vshv.tree@phylo)

#Display the tree
plot(vshv.tree@phylo)

#And let's save the tree
png(filename="imgs/1st_tree.png")
plot(vshv.tree@phylo)
dev.off()
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_07_Visualization_Pro//1st_tree.png" align='center', width="100%">
</div>

Um. Okay, but nothing special! There is a lot of work ahead to make this tree pretty!<br>

----------------------------------------------

### **Step 3: Adding host information to the tree**

Now we will read the metadata and add the host information to the tree. To do so we will use the `%<+%` sign, which allow us to use different columns of the metadata to decorate different parts of the tree.<br>

Let's add host information to the tip points.<br>

**_Input_**

```r
#Read the metadata
meta <- read.table('data/metadata.csv', sep=',', header=T)

#Plot the tree and add host information from the metadata to the tree's tip points
tree_2 <- ggtree(vshv.tree) %<+% meta +
  geom_tippoint(aes(color=Host))

#Display the tree
tree_2
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_07_Visualization_Pro//2nd_tree.png" align='center', width="100%">
</div>

And let's save the tree.<br>

**_Input_**

```r
ggsave("imgs/2nd_tree.png", tree_2, dpi = 600)
```

----------------------------------------------

### **Step 4: Use the circular layout and add bootstrap information**

`IQ-TREE` Newick file contains a bootstrap support value which has both `SH_aLRT` and `UFboot` values. Let's show branches that has high support for both of parameters (i.e. ≥ 70).<br>

When we define a `ggtree` object using `tree_3 <- ggtree(vshv.tree) %<+% meta + geom_tippoint(aes(color=Host))`, it attaches the metadata as well as tree-branch and tips info coming from original Newick tree-file in the data placeholder, which we can use by `tree_3$data`. So let’s update this to have a new column of a variable that satisfy both bootstrap values:<br>

**_Input_**

```r
#Use the circular layout for the tree
tree_3 <- ggtree(vshv.tree, layout="circular") %<+% meta +
  geom_tippoint(aes(color=Host))

#Create a new parameter `bootstrap` with the default value of 0
tree_3$data$bootstrap <- '0'

#Assign value 1 to the tree branches that has both SH_aLRT and UFboot values higher than 70
tree_3$data[which(tree_3$data$SH_aLRT >= 70 & tree_3$data$UFboot  >= 70),]$bootstrap <- '1'

#Plot the tree and use "black" color for branches with bootstrap value = 1 (>=70) and "grey" color for branches with bootstrap value = 0 (<70)
tree_3 <- tree_3 + new_scale_color() +
  geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames(c("black", "grey"), c(T,F)), guide = "none")

#Display the tree
tree_3
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_07_Visualization_Pro//3rd_tree.png" align='center', width="100%">
</div>

And let's save the tree.<br>

**_Input_**

```r
ggsave("imgs/3rd_tree.png", tree_3, dpi = 600)
```

----------------------------------------------

### **Step 5: Add the last metadata info as the heatmaps and make the final visualization**

Create two separate dataframes for `Water` and `Year` and parse them to make further visualize much more convenient.<br>

**_Input_**

```r
#Read Water column from metadata to a separate dataframe
meta.water <- as.data.frame(meta[,'Water'])
#Change column name from `meta[,'Water']` to `Water`
colnames(meta.water) <- 'Water'
#Change row names from `1, 2, 3, ...` to strains `AU-8-95, CH-FI262BFH, ...`
rownames(meta.water) <- meta$Strain

#Read Year column from metadata to a separate dataframe
meta.year <- as.data.frame(meta[,'Year'])
##Change column name from `meta[,'Year']` to `Year`
colnames(meta.year) <- 'Year'
#Change row names from `1, 2, 3, ...` to strains `AU-8-95, CH-FI262BFH, ...`
rownames(meta.year) <- meta$Strain
```

**_Input_**

```r
#Create a new tree that uses the previous tree (circular layout + host info at the tip points)
#Add `meta.water` dataframe as the heatmap and set the width to 0.2 and offset to 0.01
#Use `viridis` "A" colormap and name legend title "Water"
final_tree <- gheatmap(tree_3, meta.water, width=0.2, offset=0.01) + 
  scale_fill_viridis_d(option="A", name="Water") +
  new_scale_fill() #define a new fill scale using after the first gheatmap, for the second gheatmap to draw upon

#Create a new tree that uses the tree created above
#Add `meta.year` dataframe as the heatmap and set the width to 0.1 and offset to 0.05
#Use `viridis` "D" colormap and name legend title "Year"
final_tree <- gheatmap(final_tree, meta.year, width=0.1, offset=0.05) + 
  scale_fill_viridis(option="D", name="Year")
```

Add the country information from the metadata to the tip labels and plot the final tree!<br>

**_Input_**

```r
#Add the country information to the tip labels, set the color to "gray40", offset to 0.03 and size to 3
final_tree <- final_tree + geom_tiplab2(aes(label=Country), color="gray40", offset=0.003, size=3)

#Display the tree
final_tree
```

**_Output_**

<div style='justify-content: center'>
<img src="/images/ngs-handbook/04_Phylogenetics/04_07_Visualization_Pro//final_tree.png" align='center', width="100%">
</div>

And let's save the tree.<br>

**_Input_**

```r
ggsave("imgs/final_tree.png", final_tree, width = 10, height = 8, dpi = 600)
```