---
icon: lucide/dna
---

# **Whole Genome and Pangenome Analyses**
>This chapter contains a manual on whole genome and pangenome analyses.

We will use the `PanACoTA` tool to perform the pangenome analysis.

## **Instruction**
![CLT](https://img.shields.io/badge/Language-Command Line Tools & Python-steelblue)<br>
![IDE](https://img.shields.io/badge/Recommended IDE-Jupyter Notebook-steelblue)

Every step of this pipeline can be done in `terminal` except **Step 6.1.**, which needs `Python`.<br>
It is recommended to use `Jupyter Notebook`. Just write `!` in the beggining of each cell to make it understand `bash` commands.

To recreate any of the steps of this manual please install:

```bash
wget https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/envs/panacota.yaml
```

```bash
conda env create -f panacota.yaml
```

And of cource do not forget to activate the envinronment!

```bash
conda activate panacota
```

----------------------------------------------

## **Introduction**

Here we will work with `PanACoTA` pipeline. This tool has many dependecies:

- For `prepare` module: `mash` (to filter genomes)
- For `annotate` module: `prokka` and/or `prodigal` (to uniformly annotate your genomes)
- For `pangenome` module: `mmseqs` (to generate pangenomes)
- For `align` module: `mafft` (to align persistent genome)
- For `tree` module: At least one of those softwares, to infer a phylogenetic tree:
    - `IQ Tree`
    - `FastTreeMP`
    - `FastME`
    - `Quicktree`

All of them but one will be installed with the `panacota.yaml` conda environment.<br>
The last one - `mash` please install yourself.<br>
For more details please visit [`PanACoTA` GitHub repository](https://github.com/gem-pasteur/PanACoTA)

## **Task**

???+ info

    We've sequenced and assembled the genome.<br>
    It's the genome of some organism, and we even know what kind of organism it is.<br>

    - We download the genomes of its close relatives from public databases.<br>
    - We annotate them all together (so that there is a homogeneous annotation).<br>
    - Build orthologous series.<br>
    - Extract the pangenome.<br>
    - Make genes alignment.<br>
    - Build a phylogenetic tree.<br>

    We'll be working with the _Shigella flexneri_ genome.<br>
    It's genetically nested within _Escherichia_, but it's a very aggressive pathogen.<br>
    The same pathogenicity factors are found in _E. marmotae_ isolated from other mammals.<br>
    So our target is _Shigella flexneri_. And for comparison, we want to download the genomes of _E. marmotae_ isolated from marmosets and analyse them further.

----------------------------------------------

## **Step 1: Download genomes from the database and check their quality**

As it was mention in the **Introduction** we will run `PanACoTA` pipeline.
It has many dependecies and in the `panacota.yaml` cona environment all but one are already installed.
The one missing is `mash`.
It is needed for this step, but it is not available in `conda` or `pip`.
It is distributed as binary. To install it on `Ubuntu` run this command in your terminal:

```
sudo apt install mash
```

For more details read the [`mash` tutorial](https://mash.readthedocs.io/en/latest/)

`Mash` measures the distances between genomes:<br>

- if the distances are any small (i.e. we have some copies).<br>
- if the distances are too large (i.e. genomes so far away that nothing can be counted on them further).<br>

The programme discards them (this is a configurable parameter).<br>

???+ info
    - `GenBank` - in general, all possible `.fasta` files uploaded to `NCBI's GenBank` database.<br>
    - `RefSeq` - an attempt to structure `GenBank`, a so-called ‘quality mark’ meaning that `NCBI` staff have analysed the data according to their own pipeline and approved its quality at a high level.<br>

Let's break down the command below in parts:<br>

- `PanACoTA` - the name of the tool we are running.<br>
- `prepare` - one of the modules of this tool.<br>
- By key `-g` we pass what genomes we are downloading - quotes are very important.<br>
- By key `-s` we choose between GenBank and RefSeq.<br>
- By key `-l` we choose the level of assembly.<br>

**_Input_**

```bash
PanACoTA prepare -g 'Escherichia marmotae' -s refseq -l complete
```

Let's see what we have after executing this command

**_Input_**

```bash
ls Escherichia_marmotae
```

**_Output_**

```
Database_init
LSTINFO-Escherichia_marmotae-filtered-0.0001_0.06.txt
PanACoTA_prepare_Escherichia_marmotae.log
PanACoTA_prepare_Escherichia_marmotae.log.details
PanACoTA_prepare_Escherichia_marmotae.log.err
assembly_summary-Escherichia_marmotae.txt
discarded-by-L90_nbcont-Escherichia_marmotae.lst
discarded-by-minhash-Escherichia_marmotae-0.0001_0.06.txt
mash_files
refseq
tmp_files
```

- `PanACoTA_prepare_Escherichia_marmotae.log` - file with log of command execution.<br>
- `PanACoTA_prepare_Escherichia_marmotae.log.details` - file with detailed log of command execution.<br>
- `PanACoTA_prepare_Escherichia_marmotae.log.err` - error log (if there were any), it is always created, but if there were no errors, the file will be empty.<br>
- `assembly_summary-Escherichia_marmotae.txt` - statistics on a specific assembly of _Escherichia marmotae_.<br>
- `discarded-by-minhash-Escherichia_marmotae-0.0001_0.06.txt` - file with what `mash` discarded (it discarded anything closer than `0.0001` and anything further than `0.06`).<br>
- `LSTINFO-Escherichia_marmotae-filtered-0.0001e_0.06.txt` - a list of all the files that `PanACoTA` fragmented.<br>

Now let's look in the `Database_init` directory.<br>

**_Input_**

```bash
ls Escherichia_marmotae/Database_init
```

**_Output_**

```
GCF_002900365.1_ASM290036v1_genomic.fna
GCF_013636045.1_ASM1363604v1_genomic.fna
GCF_013636235.1_ASM1363623v1_genomic.fna
GCF_013732895.1_ASM1373289v1_genomic.fna
GCF_013745515.1_ASM1374551v1_genomic.fna
GCF_013746655.1_ASM1374665v1_genomic.fna
GCF_022592155.1_ASM2259215v1_genomic.fna
GCF_029717905.1_ASM2971790v1_genomic.fna
GCF_029719265.1_ASM2971926v1_genomic.fna
GCF_029962465.1_ASM2996246v1_genomic.fna
GCF_037055335.1_ASM3705533v1_genomic.fna
GCF_900636405.1_41767_E01_genomic.fna
GCF_900637015.1_46514_C01_genomic.fna
```
Thus `PanACoTA` downloaded 12 complete genomes (nucleotide genomic sequences) of _Escherichia marmotae_ from `RefSeq` into the `Escherichia_marmotae/Database_init directory`.<br>

Now we will add to all this our _Shigella flexneri_.

**_Input_**

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/405/GCF_000007405.1_ASM740v1/GCF_000007405.1_ASM740v1_genomic.fna.gz -P Escherichia_marmotae/Database_init/
```

Let's unpack the archive.<br>

**_Input_**

```bash
gzip -d Escherichia_marmotae/Database_init/GCF_000007405.1_ASM740v1_genomic.fna.gz
```

The next step of `PanACoTA` takes as input a text file with a list of files containing genomes.<br>
That's why we'll create it!

**_Input_**

```bash
ls Escherichia_marmotae/Database_init/ > listFile 
```

Let's check this list

**_Input_**

```bash
cat listFile
```

**_Output_**

```
GCF_000007405.1_ASM740v1_genomic.fna
GCF_002900365.1_ASM290036v1_genomic.fna
GCF_013636045.1_ASM1363604v1_genomic.fna
GCF_013636235.1_ASM1363623v1_genomic.fna
GCF_013732895.1_ASM1373289v1_genomic.fna
GCF_013745515.1_ASM1374551v1_genomic.fna
GCF_013746655.1_ASM1374665v1_genomic.fna
GCF_022592155.1_ASM2259215v1_genomic.fna
GCF_029717905.1_ASM2971790v1_genomic.fna
GCF_029719265.1_ASM2971926v1_genomic.fna
GCF_029962465.1_ASM2996246v1_genomic.fna
GCF_037055335.1_ASM3705533v1_genomic.fna
GCF_900636405.1_41767_E01_genomic.fna
GCF_900637015.1_46514_C01_genomic.fna
```

That's right!<br>
These are the files that will be given as input in the next step of `PanACoTA` pipeline.<br>

----------------------------------------------

## **Step 2: Annotate the genomes**

Running `PanACoTA` with the `annotate` mode.<br>
This mode of `PanACoTA` is basically a shell for using the `prokka` tool.<br>
Let's also break down the command below in parts:<br>

- `PanACoTA` - the name of the tool we are running.<br>
- `annotate` - module of this tool that annotates genomes.<br>
- By key `-d` we pass the directory with the genomes.<br>
- By key `-r` we specify where to store the output files.<br>
- By key `-n` we set some four-character name for the data, which will then be used in the pipeline.<br>
- By key `-l` we specify what exactly genomes from the directory in `-d` key must be annotated.<br>
- By key `--threads` we specify how much threads this command can use (if use `0` - automatically uses all available cores).<br>

???+ tip
    It takes a long time to execute the command (it took me 19 minutes)

**_Input_**

```bash
PanACoTA annotate -d Escherichia_marmotae/Database_init/ -r Annotation -n EsMa -l listFile --threads 24
```

Let's see what we got as a result of the annotation

**_Input_**

```bash
ls Annotation
```

**_Output_**

```
Genes					QC_L90-listFile.png
LSTINFO					QC_nb-contigs-listFile.png
LSTINFO-.lst				Replicons
PanACoTA-annotate_listFile.log		discarded-.lst
PanACoTA-annotate_listFile.log.details	gff3
PanACoTA-annotate_listFile.log.err	tmp_files
Proteins
```

We now have an `Annotation` folder, and in it:<br>

- `Proteins`<br>
- `Genes`<br>
- `Replicons`<br>
- `gff3`<br>
- `LSTINFO`<br>

### **Step 2.1. `Proteins` directory**

Let's look at the `Annotation/Proteins/` directory. What is stored in it?

**_Input_**

```bash
ls Annotation/Proteins/
```

**_Output_**

```
EsMa.0824.00001.prt  EsMa.0824.00006.prt  EsMa.0824.00011.prt
EsMa.0824.00002.prt  EsMa.0824.00007.prt  EsMa.0824.00012.prt
EsMa.0824.00003.prt  EsMa.0824.00008.prt  EsMa.0824.00013.prt
EsMa.0824.00004.prt  EsMa.0824.00009.prt  EsMa.0824.00014.prt
EsMa.0824.00005.prt  EsMa.0824.00010.prt  EsMa.All.prt
```

So, we have 15 files in this folder (14 files for each organism and 1 common file).<br>

- `EsMa` is the name we set in the `-n` key.<br>
- `0824` is the launch ID - August 24.<br>
- `00001-00014` - organism sequence number (IMPORTANT: the files are not always in the exact order in which they were uploaded to listFile).<br>

Let's take a look at `Annotation/Proteins/EsMa.0824.00001.prt`.

**_Input_**

```bash
head -10 Annotation/Proteins/EsMa.0824.00001.prt
```

**_Output_**

```
>EsMa.0824.00001.0001b_00001 2463 thrA | Bifunctional aspartokinase/homoserine dehydrogenase 1 | NA | similar to AA sequence:UniProtKB:P00561 | COG:COG0460
MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDA
LPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINA
ALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIP
ADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQV
PDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRD
EDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISF
CVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAAL
ARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGAL
LEQLKRQQSWLKNKHIDLRVCGVANSKALLTSVHGLNLENWQEELAQAKEPFNLGRLIRL
```

And inside each `.prt` file, we see the amino acid sequences of all the proteins that `PanACoTA` predicted for us!<br>
Let's turn our attention to line 1:<br>

**_Input_**

```bash
head -1 Annotation/Proteins/EsMa.0824.00001.prt
```

**_Output_**

```
>EsMa.0824.00001.0001b_00001 2463 thrA | Bifunctional aspartokinase/homoserine dehydrogenase 1 | NA | similar to AA sequence:UniProtKB:P00561 | COG:COG0460
```

- `>EsMa.0824.00001.` - the name of the genome in the `PanACoTA`'s annotation.<br>
- `0001b` - replicon number.<br>
- `_00001` - frame number.<br>
- `2463` - length.<br>
- `thrA` - protein's name.<br>
- `Bifunctional aspartokinase/homoserine dehydrogenase 1` - protein's function.<br>
- `similar to AA sequence:UniProtKB:P00561` - the basis for predicting the protein and its function.<br>
- `COG:COG0460` - COG ID.<br>

### **Step 2.2. `LSTINFO-.lst` file**

Now, let's look at the `Annotation/LSTINFO-.lst` file.

**_Input_**

```bash
cat Annotation/LSTINFO-.lst
```

**_Output_**

```
gembase_name	orig_name	to_annotate	gsize	nb_conts	L90
EsMa.0824.00001	GCF_000007405.1_ASM740v1_genomic.fna	Annotation/tmp_files/GCF_000007405.1_ASM740v1_genomic.fna_prokka-split5N.fna	4599354	1	1
EsMa.0824.00002	GCF_900636405.1_41767_E01_genomic.fna	Annotation/tmp_files/GCF_900636405.1_41767_E01_genomic.fna_prokka-split5N.fna	4857140	1	1
EsMa.0824.00003	GCF_900637015.1_46514_C01_genomic.fna	Annotation/tmp_files/GCF_900637015.1_46514_C01_genomic.fna_prokka-split5N.fna	4450344	1	1
EsMa.0824.00004	GCF_013636045.1_ASM1363604v1_genomic.fna	Annotation/tmp_files/GCF_013636045.1_ASM1363604v1_genomic.fna_prokka-split5N.fna	4637518	2	1
EsMa.0824.00005	GCF_029962465.1_ASM2996246v1_genomic.fna	Annotation/tmp_files/GCF_029962465.1_ASM2996246v1_genomic.fna_prokka-split5N.fna	4669757	2	1
EsMa.0824.00006	GCF_002900365.1_ASM290036v1_genomic.fna	Annotation/tmp_files/GCF_002900365.1_ASM290036v1_genomic.fna_prokka-split5N.fna	4896291	3	1
EsMa.0824.00007	GCF_013636235.1_ASM1363623v1_genomic.fna	Annotation/tmp_files/GCF_013636235.1_ASM1363623v1_genomic.fna_prokka-split5N.fna	4981845	3	1
EsMa.0824.00008	GCF_013745515.1_ASM1374551v1_genomic.fna	Annotation/tmp_files/GCF_013745515.1_ASM1374551v1_genomic.fna_prokka-split5N.fna	4887958	3	1
EsMa.0824.00009	GCF_022592155.1_ASM2259215v1_genomic.fna	Annotation/tmp_files/GCF_022592155.1_ASM2259215v1_genomic.fna_prokka-split5N.fna	5200167	3	1
EsMa.0824.00010	GCF_013746655.1_ASM1374665v1_genomic.fna	Annotation/tmp_files/GCF_013746655.1_ASM1374665v1_genomic.fna_prokka-split5N.fna	5094503	4	1
EsMa.0824.00011	GCF_037055335.1_ASM3705533v1_genomic.fna	Annotation/tmp_files/GCF_037055335.1_ASM3705533v1_genomic.fna_prokka-split5N.fna	5022908	4	1
EsMa.0824.00012	GCF_029719265.1_ASM2971926v1_genomic.fna	Annotation/tmp_files/GCF_029719265.1_ASM2971926v1_genomic.fna_prokka-split5N.fna	5075061	5	1
EsMa.0824.00013	GCF_013732895.1_ASM1373289v1_genomic.fna	Annotation/tmp_files/GCF_013732895.1_ASM1373289v1_genomic.fna_prokka-split5N.fna	5083937	6	1
EsMa.0824.00014	GCF_029717905.1_ASM2971790v1_genomic.fna	Annotation/tmp_files/GCF_029717905.1_ASM2971790v1_genomic.fna_prokka-split5N.fna	4780041	7	1
```

- `EsMa.0824.00001` - the name that `PanACoTA` renamed the file to (by `-n` key).<br>
- `GCF_000007405.1_ASM740v1_genomic.fna` - the name from which `PanACoTA` was renaming.<br>
- `Annotation/tmp_files/GCF_000007405.1_ASM740v1_genomic.fna_prokka-split5N.fna`- a folder containing what `PanACoTA` eventually annotated.<br>
- `4599354` - genome size.<br>

### **Step 2.3. `Genes` directory**

Let's take a look at the contents of the `Genes` directory.<br>

**_Input_**

```bash
ls Annotation/Genes/
```

**_Output_**

```
EsMa.0824.00001.gen  EsMa.0824.00006.gen  EsMa.0824.00011.gen
EsMa.0824.00002.gen  EsMa.0824.00007.gen  EsMa.0824.00012.gen
EsMa.0824.00003.gen  EsMa.0824.00008.gen  EsMa.0824.00013.gen
EsMa.0824.00004.gen  EsMa.0824.00009.gen  EsMa.0824.00014.gen
EsMa.0824.00005.gen  EsMa.0824.00010.gen
```

So we have 14 files in this folder (14 files for each organism).<br>

**_Input_**

```bash
head -10 Annotation/Genes/EsMa.0824.00001.gen
```

**_Output_**

```
>EsMa.0824.00001.0001b_00001 2463 thrA | Bifunctional aspartokinase/homoserine dehydrogenase 1 | NA | similar to AA sequence:UniProtKB:P00561 | COG:COG0460
ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTT
GCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCC
GCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCT
TTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCC
GCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAA
ATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCT
GCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTGTTAGAAGCG
CGTGGTCACAACGTTACCGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTAC
CTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCG
```

So here we have nucleotide sequences! The names of genes and proteins are the same!<br>

### **Step 2.4. `Replicons` directory**

We also have a directory called `Replicons`, but it's a boring directory that just holds actually the original nucleotide sequences, but just under code numbers specified in the `-n` key.<br>

**_Input_**

```bash
ls Annotation/Replicons/
```

**_Output_**

```
EsMa.0824.00001.fna  EsMa.0824.00006.fna  EsMa.0824.00011.fna
EsMa.0824.00002.fna  EsMa.0824.00007.fna  EsMa.0824.00012.fna
EsMa.0824.00003.fna  EsMa.0824.00008.fna  EsMa.0824.00013.fna
EsMa.0824.00004.fna  EsMa.0824.00009.fna  EsMa.0824.00014.fna
EsMa.0824.00005.fna  EsMa.0824.00010.fna
```

**_Input_**

```bash
head Annotation/Replicons/EsMa.0824.00001.fna
```

**_Output_**

```
>EsMa.0824.00001.0001 4599354
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC
TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG
TCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTAC
ACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT
AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGG
CTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTA
CATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCA
GGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGG
CGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAAC
```

### **Step 2.5. `LSTINFO` directory**

And we have a directory called `LSTINFO`.<br>

**_Input_**

```bash
ls Annotation/LSTINFO/
```

**_Output_**

```
EsMa.0824.00001.lst  EsMa.0824.00006.lst  EsMa.0824.00011.lst
EsMa.0824.00002.lst  EsMa.0824.00007.lst  EsMa.0824.00012.lst
EsMa.0824.00003.lst  EsMa.0824.00008.lst  EsMa.0824.00013.lst
EsMa.0824.00004.lst  EsMa.0824.00009.lst  EsMa.0824.00014.lst
EsMa.0824.00005.lst  EsMa.0824.00010.lst
```

**_Input_**

```bash
head -10 Annotation/LSTINFO/EsMa.0824.00001.lst
```

**_Output_**

```
336	2798	D	CDS	EsMa.0824.00001.0001b_00001	thrA	| Bifunctional aspartokinase/homoserine dehydrogenase 1 | NA | similar to AA sequence:UniProtKB:P00561 | COG:COG0460
2800	3732	D	CDS	EsMa.0824.00001.0001i_00002	thrB	| Homoserine kinase | 2.7.1.39 | similar to AA sequence:UniProtKB:P00547 | COG:COG0083
3733	5019	D	CDS	EsMa.0824.00001.0001i_00003	thrC	| Threonine synthase | 4.2.3.1 | similar to AA sequence:UniProtKB:P00934 | COG:COG0498
5233	5529	D	CDS	EsMa.0824.00001.0001i_00004	NA	| hypothetical protein | NA | NA | NA
5682	6458	C	CDS	EsMa.0824.00001.0001i_00005	yaaA	| Peroxide stress resistance protein YaaA | NA | similar to AA sequence:UniProtKB:P0A8I3 | COG:COG3022
6528	6896	C	CDS	EsMa.0824.00001.0001i_00006	alsT_1	| Amino-acid carrier protein AlsT | NA | similar to AA sequence:UniProtKB:Q45068 | COG:COG1115
6918	7958	C	CDS	EsMa.0824.00001.0001i_00007	alsT_2	| Amino-acid carrier protein AlsT | NA | similar to AA sequence:UniProtKB:Q45068 | COG:COG1115
8237	9190	D	CDS	EsMa.0824.00001.0001i_00008	talB	| Transaldolase B | 2.2.1.2 | similar to AA sequence:UniProtKB:P0A870 | COG:COG0176
9305	9892	D	CDS	EsMa.0824.00001.0001i_00009	mog	| Molybdopterin adenylyltransferase | 2.7.7.75 | similar to AA sequence:UniProtKB:P0AF03 | COG:COG0521
9927	10493	C	CDS	EsMa.0824.00001.0001i_00010	satP	| Succinate-acetate/proton symporter SatP | NA | similar to AA sequence:UniProtKB:P0AC98 | COG:COG1584
```

This is a tabular format of all annotated features in the organism!<br>
This is especially convenient if we work with genomes using programming languages (e.g. `Python` and `pandas`).<br>

### **Step 2.6. `tmp_files` directory**

**_Input_**

```bash
ls -U Annotation/tmp_files/ | head -10
```
**_Output_**

```
GCF_000007405.1_ASM740v1_genomic.fna_prokka-split5N.fna
GCF_000007405.1_ASM740v1_genomic.fna_prokka-split5N.fna-prokka.log
GCF_000007405.1_ASM740v1_genomic.fna_prokka-split5N.fna-prokkaRes
GCF_002900365.1_ASM290036v1_genomic.fna_prokka-split5N.fna
GCF_002900365.1_ASM290036v1_genomic.fna_prokka-split5N.fna-prokka.log
GCF_002900365.1_ASM290036v1_genomic.fna_prokka-split5N.fna-prokkaRes
GCF_013636045.1_ASM1363604v1_genomic.fna_prokka-split5N.fna
GCF_013636045.1_ASM1363604v1_genomic.fna_prokka-split5N.fna-prokka.log
GCF_013636045.1_ASM1363604v1_genomic.fna_prokka-split5N.fna-prokkaRes
GCF_013636235.1_ASM1363623v1_genomic.fna_prokka-split5N.fna
```

The `tmp_files` directory contains logs on temporary files.<br>
What is it for?<br>
If, let's say, something that should be there is not found in the genome under study, the first thing to do is to study these files - most likely `prokka` crashed with an error and you can read the details in these temporary files.<br>

----------------------------------------------

## **Step 3: Build orthological series**

This command builds the pangenome.<br>
Let's also break down the command below in parts:<br>

- `PanACoTA` - the name of the tool we are running.<br>
- `pangenome` - module of this tool that builds the pangenome.<br>
- By key `-l` we pass the `LSTINFO-.lst` file from the **Step 2.2.**.<br>
- By key `-n` we pass the four-character name for the data, which was set previously.<br>
- By key `-d` we pass the directory with the proteins (since `PanACoTA` builds the pangenome by protein sequence similarity, it only needs the protein directory).<br>
- By key `-o` we specify where to store the output files.<br>
- By the key `-i` we establish the level of similarity (identity) of the sequences.<br>

**_Input_**

```bash
PanACoTA pangenome -l Annotation/LSTINFO-.lst -n EsMa -d Annotation/Proteins/ -o Pangenome -i 0.8 
```

Let's take a look at the log file.<br>

**_Input_**

```bash
cat Pangenome/PanACoTA-pangenome_EsMa.log
```

**_Output_**

```
[2024-08-25 18:54:07] :: INFO :: PanACoTA version 1.4.0
[2024-08-25 18:54:07] :: INFO :: Command used
 	 > PanACoTA pangenome -l Annotation/LSTINFO-.lst -n EsMa -d Annotation/Proteins/ -o Pangenome -i 0.8
[2024-08-25 18:54:07] :: INFO :: Building bank with all proteins to Annotation/Proteins/EsMa.All.prt
[2024-08-25 18:54:08] :: INFO :: Will run MMseqs2 with:
	- minimum sequence identity = 80.0%
	- cluster mode 1
[2024-08-25 18:54:08] :: INFO :: Creating database
[2024-08-25 18:54:11] :: INFO :: Clustering proteins...
[2024-08-25 18:54:36] :: INFO :: Converting mmseqs results to pangenome file
[2024-08-25 18:54:36] :: INFO :: Pangenome has 9694 families.
[2024-08-25 18:54:36] :: INFO :: Retrieving information from pan families
[2024-08-25 18:54:36] :: INFO :: Generating qualitative and quantitative matrix, and summary file
[2024-08-25 18:54:37] :: INFO :: DONE
```

???+ info
    `[2024-08-25 18:54:36] :: INFO :: Pangenome has 9694 families.`<br>

That is, in this command, `PanACoTA` counted a total of `9694` protein families.

In this step, we took the annotated genomes, we took the proteins from them.<br>
And then the orthology table is constructed with them.<br>
The outcome of this step will be lists of genes grouped together.<br>
As output, we're going to see a bunch of logs, and one single table.

**_Input_**

```bash
ls Pangenome/
```

**_Output_**

```
PanACoTA-pangenome_EsMa.log
PanACoTA-pangenome_EsMa.log.details
PanACoTA-pangenome_EsMa.log.err
PanGenome-EsMa.All.prt-clust-0.8-mode1.lst
PanGenome-EsMa.All.prt-clust-0.8-mode1.lst.bin
PanGenome-EsMa.All.prt-clust-0.8-mode1.lst.quali.txt
PanGenome-EsMa.All.prt-clust-0.8-mode1.lst.quanti.txt
PanGenome-EsMa.All.prt-clust-0.8-mode1.lst.summary.txt
mmseq_EsMa.All.prt_0.8-mode1.log
tmp_EsMa.All.prt_0.8-mode1
```

The only file of interest is `PanGenome-EsMa.All.prt-clust-0.8-mode1.lst`.<br>
Let's take a look at it.

**_Input_**

```bash
head -10 Pangenome/PanGenome-EsMa.All.prt-clust-0.8-mode1.lst
```

**_Output_**

```
1 EsMa.0824.00001.0001i_04241 EsMa.0824.00002.0001i_00106 EsMa.0824.00003.0001i_00093 EsMa.0824.00004.0001i_00095 EsMa.0824.00005.0001i_00094 EsMa.0824.00006.0001i_03340 EsMa.0824.00007.0001i_02672 EsMa.0824.00008.0001i_00094 EsMa.0824.00009.0001i_00094 EsMa.0824.00010.0001i_04572 EsMa.0824.00011.0001i_02424 EsMa.0824.00012.0001i_02999 EsMa.0824.00013.0001i_02846 EsMa.0824.00014.0001i_00094
2 EsMa.0824.00001.0001i_03873 EsMa.0824.00002.0001i_00130 EsMa.0824.00003.0001i_00117 EsMa.0824.00004.0001i_00119 EsMa.0824.00005.0001i_00118 EsMa.0824.00006.0001i_03364 EsMa.0824.00007.0001i_02648 EsMa.0824.00008.0001i_00118 EsMa.0824.00009.0001i_00126 EsMa.0824.00010.0001i_04604 EsMa.0824.00011.0001i_02400 EsMa.0824.00012.0001i_03032 EsMa.0824.00013.0001i_02821 EsMa.0824.00014.0001i_00126
3 EsMa.0824.00008.0002i_04583 EsMa.0824.00010.0002i_04636 EsMa.0824.00012.0002i_04658
4 EsMa.0824.00008.0002i_04649 EsMa.0824.00010.0002i_04668 EsMa.0824.00012.0002i_04736
5 EsMa.0824.00007.0003i_04731 EsMa.0824.00008.0002i_04616 EsMa.0824.00010.0002i_04700 EsMa.0824.00012.0002i_04701 EsMa.0824.00012.0003i_04789
6 EsMa.0824.00010.0003i_04833
7 EsMa.0824.00001.0001i_02181 EsMa.0824.00002.0001i_02692 EsMa.0824.00003.0001i_02407 EsMa.0824.00004.0001i_02595 EsMa.0824.00005.0001i_02480 EsMa.0824.00006.0001i_00413 EsMa.0824.00007.0001i_00212 EsMa.0824.00008.0001i_02711 EsMa.0824.00009.0001i_02866 EsMa.0824.00010.0001i_02477 EsMa.0824.00011.0001i_00009 EsMa.0824.00011.0001i_04359 EsMa.0824.00012.0001i_01000 EsMa.0824.00013.0001i_01030 EsMa.0824.00014.0001i_02416
8 EsMa.0824.00005.0001i_02448 EsMa.0824.00011.0001i_00041 EsMa.0824.00011.0001i_04391 EsMa.0824.00011.0001i_04451
9 EsMa.0824.00002.0001i_03089 EsMa.0824.00003.0001i_02788 EsMa.0824.00005.0001i_02414 EsMa.0824.00008.0001i_02644 EsMa.0824.00008.0001i_02862 EsMa.0824.00011.0001i_00075 EsMa.0824.00011.0001i_04425 EsMa.0824.00011.0001i_04485 EsMa.0824.00014.0001i_02479
10 EsMa.0824.00001.0001i_02131 EsMa.0824.00002.0001i_02651 EsMa.0824.00003.0001i_02366 EsMa.0824.00004.0001i_02431 EsMa.0824.00005.0001i_02379 EsMa.0824.00006.0001i_00454 EsMa.0824.00007.0001i_00253 EsMa.0824.00008.0001i_02607 EsMa.0824.00009.0001i_02686 EsMa.0824.00010.0001i_02436 EsMa.0824.00011.0001i_00110 EsMa.0824.00012.0001i_00951 EsMa.0824.00013.0001i_00978 EsMa.0824.00014.0001i_02375
```

In this file, each line corresponds to one orthologous row.<br>
The identifiers of all the genes that fall into this row are written in the line with a separator.<br>
For example, in line 6 we have one gene - it means that it is a singleton (it is so one, and no one has anything similar).<br>
Long rows are orthologous groups with many genes in them.<br>

Further after this step, from the orthologous series, we identify the `core` - that is, the genes encoding proteins that occur in a single copy in each genome.<br>

----------------------------------------------

## **Step 4: Extract pangenome**

Let's also break down the command below in parts:<br>

- `PanACoTA` - the name of the tool we are running.<br>
- `corepers` - module of this tool that identifies and extracts `core`.<br>
- By key `-p` we pass the `PanGenome-EsMa.All.prt-clust-0.8-mode1.lst` file.<br>
- By key `-o` we specify where to store the output files.<br>

**_Input_**

```bash
PanACoTA corepers -p Pangenome/PanGenome-EsMa.All.prt-clust-0.8-mode1.lst -o Coregenome
```

Let's look at the log file.<br>

**_Input_**

```bash
cat Coregenome/PanACoTA-corepers.log
```

**_Output_**

```
[2024-08-25 18:54:45] :: INFO :: PanACoTA version 1.4.0
[2024-08-25 18:54:45] :: INFO :: Command used
 	 > PanACoTA corepers -p Pangenome/PanGenome-EsMa.All.prt-clust-0.8-mode1.lst -o Coregenome
[2024-08-25 18:54:45] :: INFO :: Will generate a CoreGenome.
[2024-08-25 18:54:45] :: INFO :: Retrieving info from binary file
[2024-08-25 18:54:45] :: INFO :: Generating Persistent genome of a dataset containing 14 genomes
[2024-08-25 18:54:45] :: INFO :: The core genome contains 2709 families, each one having exactly 14 members, from the 14 different genomes.
[2024-08-25 18:54:45] :: INFO :: Persistent genome step done.
```

???+ info
    `[2024-08-25 18:54:45] :: INFO :: The core genome contains 2709 families, each one having exactly 14 members, from the 14 different genomes.`<br>

That is, out of `9694` protein families, `2709` are single-copy and universal families.<br>

Now let's take a look at the output.<br>

**_Input_**

```bash
ls Coregenome/
```

**_Output_**

```
PanACoTA-corepers.log
PanACoTA-corepers.log.err
PersGenome_PanGenome-EsMa.All.prt-clust-0.8-mode1.lst-all_1.lst
```

We only have logs and one table - `PersGenome_PanGenome-EsMa.All.prt-clust-0.8-mode1.lst-all_1.lst`.<br>
Let's take a look at it.<br>

**_Input_**

```bash
head -10 Coregenome/PersGenome_PanGenome-EsMa.All.prt-clust-0.8-mode1.lst-all_1.lst
```

**_Output_**

```
1 EsMa.0824.00001.0001i_04241 EsMa.0824.00002.0001i_00106 EsMa.0824.00003.0001i_00093 EsMa.0824.00004.0001i_00095 EsMa.0824.00005.0001i_00094 EsMa.0824.00006.0001i_03340 EsMa.0824.00007.0001i_02672 EsMa.0824.00008.0001i_00094 EsMa.0824.00009.0001i_00094 EsMa.0824.00010.0001i_04572 EsMa.0824.00011.0001i_02424 EsMa.0824.00012.0001i_02999 EsMa.0824.00013.0001i_02846 EsMa.0824.00014.0001i_00094
2 EsMa.0824.00001.0001i_03873 EsMa.0824.00002.0001i_00130 EsMa.0824.00003.0001i_00117 EsMa.0824.00004.0001i_00119 EsMa.0824.00005.0001i_00118 EsMa.0824.00006.0001i_03364 EsMa.0824.00007.0001i_02648 EsMa.0824.00008.0001i_00118 EsMa.0824.00009.0001i_00126 EsMa.0824.00010.0001i_04604 EsMa.0824.00011.0001i_02400 EsMa.0824.00012.0001i_03032 EsMa.0824.00013.0001i_02821 EsMa.0824.00014.0001i_00126
10 EsMa.0824.00001.0001i_02131 EsMa.0824.00002.0001i_02651 EsMa.0824.00003.0001i_02366 EsMa.0824.00004.0001i_02431 EsMa.0824.00005.0001i_02379 EsMa.0824.00006.0001i_00454 EsMa.0824.00007.0001i_00253 EsMa.0824.00008.0001i_02607 EsMa.0824.00009.0001i_02686 EsMa.0824.00010.0001i_02436 EsMa.0824.00011.0001i_00110 EsMa.0824.00012.0001i_00951 EsMa.0824.00013.0001i_00978 EsMa.0824.00014.0001i_02375
11 EsMa.0824.00001.0001i_02083 EsMa.0824.00002.0001i_02617 EsMa.0824.00003.0001i_02332 EsMa.0824.00004.0001i_02399 EsMa.0824.00005.0001i_02347 EsMa.0824.00006.0001i_00487 EsMa.0824.00007.0001i_00287 EsMa.0824.00008.0001i_02575 EsMa.0824.00009.0001i_02653 EsMa.0824.00010.0001i_02404 EsMa.0824.00011.0001i_00142 EsMa.0824.00012.0001i_00918 EsMa.0824.00013.0001i_00944 EsMa.0824.00014.0001i_02342
12 EsMa.0824.00001.0001i_01323 EsMa.0824.00002.0001i_02588 EsMa.0824.00003.0001i_02303 EsMa.0824.00004.0001i_02370 EsMa.0824.00005.0001i_02318 EsMa.0824.00006.0001i_01160 EsMa.0824.00007.0001i_00316 EsMa.0824.00008.0001i_02546 EsMa.0824.00009.0001i_02534 EsMa.0824.00010.0001i_02375 EsMa.0824.00011.0001i_00174 EsMa.0824.00012.0001i_00889 EsMa.0824.00013.0001i_00915 EsMa.0824.00014.0001i_02313
14 EsMa.0824.00001.0001i_01389 EsMa.0824.00002.0001i_02519 EsMa.0824.00003.0001i_02241 EsMa.0824.00004.0001i_02307 EsMa.0824.00005.0001i_02256 EsMa.0824.00006.0001i_01096 EsMa.0824.00007.0001i_00379 EsMa.0824.00008.0001i_02424 EsMa.0824.00009.0001i_02389 EsMa.0824.00010.0001i_02244 EsMa.0824.00011.0001i_00240 EsMa.0824.00012.0001i_00759 EsMa.0824.00013.0001i_00851 EsMa.0824.00014.0001i_02250
15 EsMa.0824.00001.0001i_01420 EsMa.0824.00002.0001i_02486 EsMa.0824.00003.0001i_02209 EsMa.0824.00004.0001i_02275 EsMa.0824.00005.0001i_02224 EsMa.0824.00006.0001i_01062 EsMa.0824.00007.0001i_00411 EsMa.0824.00008.0001i_02391 EsMa.0824.00009.0001i_02356 EsMa.0824.00010.0001i_02211 EsMa.0824.00011.0001i_00272 EsMa.0824.00012.0001i_00727 EsMa.0824.00013.0001i_00819 EsMa.0824.00014.0001i_02217
16 EsMa.0824.00001.0001i_01485 EsMa.0824.00002.0001i_02455 EsMa.0824.00003.0001i_02177 EsMa.0824.00004.0001i_02244 EsMa.0824.00005.0001i_02192 EsMa.0824.00006.0001i_01030 EsMa.0824.00007.0001i_00442 EsMa.0824.00008.0001i_02360 EsMa.0824.00009.0001i_02325 EsMa.0824.00010.0001i_02180 EsMa.0824.00011.0001i_00304 EsMa.0824.00012.0001i_00696 EsMa.0824.00013.0001i_00787 EsMa.0824.00014.0001i_02174
20 EsMa.0824.00001.0001i_01791 EsMa.0824.00002.0001i_02262 EsMa.0824.00003.0001i_01989 EsMa.0824.00004.0001i_02098 EsMa.0824.00005.0001i_01999 EsMa.0824.00006.0001i_00822 EsMa.0824.00007.0001i_00574 EsMa.0824.00008.0001i_02152 EsMa.0824.00009.0001i_02124 EsMa.0824.00010.0001i_01983 EsMa.0824.00011.0001i_00496 EsMa.0824.00012.0001i_00561 EsMa.0824.00013.0001i_00575 EsMa.0824.00014.0001i_02021
22 EsMa.0824.00001.0001i_01879 EsMa.0824.00002.0001i_02202 EsMa.0824.00003.0001i_01929 EsMa.0824.00004.0001i_02038 EsMa.0824.00005.0001i_01936 EsMa.0824.00006.0001i_00757 EsMa.0824.00007.0001i_00637 EsMa.0824.00008.0001i_02089 EsMa.0824.00009.0001i_02064 EsMa.0824.00010.0001i_01923 EsMa.0824.00011.0001i_00560 EsMa.0824.00012.0001i_00501 EsMa.0824.00013.0001i_00512 EsMa.0824.00014.0001i_01961
```

So here we are left with only those genes, those orthologues, which are in a single copy in all the organisms we have analysed.<br>

At this step, it is possible to build more than just coregens.<br>
The `corepers` module has many more arguments, such as:<br>

- `-t` <tol>: % (between 0 and 1) of the persistent genome: a family is considered as persistent if it contains exactly one member in at least tol% of the genomes, and is absent in all other genomes. Default value for t is 1, meaning that all genomes must have a unique member. This corresponds to the coregenome (so no need to put this option if you want a coregenome). More relaxed definitions of a persistent genome can be used by using -X or -M options (see below).<br>
- `-X`: add this option if you want to relax a little the definition of the persistent genome, to get a mixed persistent genome. With -X option, a family is considered as persistent if at least tol% (tol defined by -t <tol> parameter, see above) of the genomes have exactly one member in the family, but the other genomes can have either 0, either several members in the family. This is useful to add the families where, in some genomes, 1 protein has been split in several parts, because of sequencing or assembly error(s).<br>
- `-M`: not compatible with -X. With this option, you get the multi persistent genome. It includes the strict and mixed persistent, but is even wider: the only condition for a family to be persistent is that it must have at least one member in at least tol% (tol still defined by -t <tol> parameter) of the genomes (independent of the copy number).<br>

----------------------------------------------

## **Step 5: Aligning the sequence of core gene sequences**

???+ info
    For each of our objects, we have an amino acid sequence and a nucleotide sequence. That is, we can do alignments of both.<br>
    For biological analysis, the protein sequence is preferred because:<br>

    1. Larger alphabet = more letter variants - 20 vs. 4 for nucleotide sequences.<br>
    2. Amino acids have a more conservative structure (nucleotides are more variable, there are nonsynonymous substitutions, and amino acids are under stronger selection and it is easier for us to match conservative parts).<br>

    However, when aligning amino acid sequences, we find it difficult to ‘go back’ and account for nucleotide substitutions.<br>
    But there is a trick - we take the nucleotide sequences, translate them into amino acids and sort of align the amino acids, preserving the original triplets of nucleotides that were underneath those amino acids.<br>

    Then:<br>

    - We can account for honest correct reading frame shifts - like a single nucleotide dropout in a nucleotide alignment (i.e. we allow for a situation that is impossible in an amino acid alignment).<br>
    - And at the same time we have information about the amino acid structure of the protein - a combination of two levels of information.<br>

Let's also break down the command below in parts:<br>

- `PanACoTA` - the name of the tool we are running.<br>
- `align` - module of this tool that aligns sequences with `mafft`.<br>
- By key `-c` we pass the `PersGenome_PanGenome-EsMa.All.prt-clust-0.8-mode1.lst-all_1.lst` file.<br>
- By key `-l` we pass the `LSTINFO-.lst` file from the **Step 2.2.**.<br>
- By key `-n` we pass the four-character name for the data, which was set previously.<br>
- By key `-d` we pass the path to the output directory from **Step 2.**.<br>
- By key `-o` we specify where to store the output files.<br>

**_Input_**

```bash
PanACoTA align -c Coregenome/PersGenome_PanGenome-EsMa.All.prt-clust-0.8-mode1.lst-all_1.lst -l Annotation/LSTINFO-.lst -n EsMa -d Annotation/ -o Alignment
```

Let's go into the directory and look at the results.<br>

**_Input_**

```
ls Alignment
```

**_Output_**

```
Align-EsMa  PanACoTA-align_EsMa.log	     PanACoTA-align_EsMa.log.err
List-EsMa   PanACoTA-align_EsMa.log.details  Phylo-EsMa
```

We see a lot of log files.<br>
But now we are interested in the contents of the `Align-EsMa` directory.<br>

#### **Step 5.1. `Align-EsMa` directory**

**_Input_**

```
ls -U Alignment/Align-EsMa | head -20
```

**_Output_**

```
EsMa-complete.nucl.cat.aln
EsMa-current.1.gen
EsMa-current.1.miss.lst
EsMa-current.1.prt
EsMa-current.10.gen
EsMa-current.10.miss.lst
EsMa-current.10.prt
EsMa-current.1000.gen
EsMa-current.1000.miss.lst
EsMa-current.1000.prt
EsMa-current.1007.gen
EsMa-current.1007.miss.lst
EsMa-current.1007.prt
EsMa-current.1009.gen
EsMa-current.1009.miss.lst
EsMa-current.1009.prt
EsMa-current.1011.gen
EsMa-current.1011.miss.lst
EsMa-current.1011.prt
EsMa-current.1012.gen
ls: write error: Broken pipe
```

There are a lot of files in this folder (more than 1000), so we displayed only 20 of them on the screen.<br>
But now each separate file is not one organism (as it was in the previous steps), but one orthologous series (and it contains sequences from different organisms that were glued together).<br>

- `.prt` files correspond to protein, i.e. there will be amino acid sequences inside.<br>
- `.gen` files correspond to nucleotides, so there will be nucleotide sequences inside.<br>

### **Step 5.2. `Phylo-EsMa` directory**

And now let's have a look at the contents of the `Phylo-EsMa` directory.<br>

**_Input_**

```
ls Alignment/Phylo-EsMa/
```

**_Output_**

```
EsMa.nucl.grp.aln
```

`EsMa.nucl.grp.aln` is the alignment itselft, the main result of running the `PanACoTA`'s pipeline `align` module.<br>
It is the concatenate of the core genes alignment.<br>
It is from this alignment that we will build the phylogenetic tree in the next step.<br>

----------------------------------------------

## **Step 6: Build a phylogenomic tree**

Let's also break down the command below in parts:<br>

- `PanACoTA` - the name of the tool we are running.<br>
- `tree` - module of this tool for phylogenetic tree construction.<br>
- By key `-a` we pass the alignment file (`EsMa.nucl.grp.aln`).<br>
- By key `-s` we specify which tree building tool we want to use.<br>
- By key `-o` we specify where to store the output files.<br>

**_Input_**

```bash
PanACoTA tree -a Alignment/Phylo-EsMa/EsMa.nucl.grp.aln -s iqtree2 -o Tree
```

Let's look at the output!<br>

**_Input_**

```bash
ls Tree
```

**_Output_**

```
EsMa.nucl.grp.aln.iqtree_tree.bionj   EsMa.nucl.grp.aln.iqtree_tree.treefile
EsMa.nucl.grp.aln.iqtree_tree.ckp.gz  PanACoTA-tree-iqtree2.log
EsMa.nucl.grp.aln.iqtree_tree.iqtree  PanACoTA-tree-iqtree2.log.details
EsMa.nucl.grp.aln.iqtree_tree.log     PanACoTA-tree-iqtree2.log.err
EsMa.nucl.grp.aln.iqtree_tree.mldist
```

- `PanACoTA-tree-iqtree2.log`.<br>
- `PanACoTA-tree-iqtree2.log.details`.<br>
- `PanACoTA-tree-iqtree2.log.err`.<br>

It's all `PanACoTA`'s own logs

- `EsMa.nucl.grp.aln.iqtree_tree.bionj` - this file contains the initial tree representation (tree draft).<br>
- `EsMa.nucl.grp.aln.iqtree_tree.treefile` - the very tree file we are interested in.<br>

Let's take a look at it!<br>

**_Input_**

```bash
cat Tree/EsMa.nucl.grp.aln.iqtree_tree.treefile
```

**_Output_**

```
(EsMa.0824.00001:0.0899533237,((EsMa.0824.00002:0.0014558718,(((EsMa.0824.00003:0.0012307785,EsMa.0824.00012:0.0012870353):0.0002683664,(EsMa.0824.00004:0.0013314298,(EsMa.0824.00006:0.0026908476,((EsMa.0824.00009:0.0000384090,EsMa.0824.00010:0.0000274042):0.0015596552,EsMa.0824.00013:0.0014403274):0.0002421241):0.0003947757):0.0002574774):0.0000880096,(EsMa.0824.00007:0.0012606069,EsMa.0824.00008:0.0011876549):0.0002556339):0.0001349309):0.0002322047,EsMa.0824.00014:0.0012843745):0.0012230243,(EsMa.0824.00005:0.0000163401,EsMa.0824.00011:0.0001690095):0.0009982711);
```

This is a tree in `Newick` format.<br>

### **Step 6.1. Tree visualization**

For more details please visit [Phylogenetics handbook](https://github.com/ilypopv/NGS-Handbook/tree/main/04_Phylogenetics).<br>
Here we will not use `iTOL` or `FigTree` but we will visualize the tree in the fastest way possible - using pseudo-graphics.<br>

Import `Phylo` module from `Biopython`.<br>

**_Input_**

```python
from Bio import Phylo
```

Read the `EsMa.nucl.grp.aln.iqtree_tree.treefile` file to the `tree` variable.<br>

**_Input_**

```python
tree = Phylo.read("Tree/EsMa.nucl.grp.aln.iqtree_tree.treefile", "newick")
```

Visualize the tree!<br>

**_Input_**

```python
Phylo.draw_ascii(tree)
```

**_Output_**

```
  __________________________________________________________ EsMa.0824.00001
 |
 , EsMa.0824.00002
 |
 |, EsMa.0824.00003
 ||
 || EsMa.0824.00012
 ||
 |, EsMa.0824.00004
 ||
 ||_ EsMa.0824.00006
 ||
_||, EsMa.0824.00009
 |,|
 ||| EsMa.0824.00010
 ||
 || EsMa.0824.00013
 ||
 |, EsMa.0824.00007
 ||
 || EsMa.0824.00008
 |
 | EsMa.0824.00014
 |
 , EsMa.0824.00005
 |
 | EsMa.0824.00011
```

 That's right! The longest branch is _Shigella flexneri_.<br>