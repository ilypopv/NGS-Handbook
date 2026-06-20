---
icon: lucide/briefcase-medical
---

# **Genomic Variation Analysis**
>This chapter contains a manual on variant calling pipeline. Details are described here only on how to conduct such a study from the time the data are downloaded to the time the genetic variant annotation files are received. The biological interpretation of the results is not given here.

## **Instruction**
![CLT](https://img.shields.io/badge/Language-Command Line Tools-steelblue)<br>
![IDE](https://img.shields.io/badge/Recommended IDE-Jupyter Notebook-steelblue)

You can run commands below in your `terminal`.<br>
But if you want to write a convenient to read laboratory journal you can use `Jupyter Notebook`.<br>
In that case write `!` in the beggining of each cell to make it understand `bash` commands.

To recreate any of the steps of this manual please install:

```bash
wget https://github.com/ilypopv/NGS-Handbook/raw/refs/heads/main/envs/varcall.yaml
```

```bash
conda env create -f varcall.yaml
```

And of cource do not forget to activate the envinronment!

```bash
conda activate varcall
```

----------------------------------------------

### **Step 1: Downloading the data**

First we create all the directories to store the data

**_Input_**

```bash
mkdir data/
mkdir data/reference/
```

**Download reference file from NCBI in FASTA format**

**_Input_**

```bash
efetch -db nuccore -id U00096.3 -format fasta > data/reference/EcoliK12MG1655.fa
```

**Download reference file from NCBI in GB format**

**_Input_**

```bash
efetch -db nuccore -id U00096.3 -format gb > data/reference/EcoliK12MG1655.gb
```

**Download sample FASTQ files from SRA via `sra-tools`**

**_Input_**

```bash
fastq-dump -v --split-3 --gzip SRR17909485 -O data
```

**_Output_**

```
Preference setting is: Prefer SRA Normalized Format files with full base quality scores if available.
SRR17909485 is an SRA Normalized Format file with full base quality scores.
Read 1251776 spots for SRR17909485
Written 1251776 spots for SRR17909485
```

----------------------------------------------

### **Step 2: Check the reads quality**

**_Input_**

```bash
mkdir data/fastqc_results/
```

**_Input_**

```bash
fastqc data/SRR17909485_*.fastq.gz -o data/fastqc_results
```

**_Output_**

Results are saved to:<br>
- `SRR17909485_1_fastqc.html`<br>
- `SRR17909485_2_fastqc.html`<br>
Please open them in any web browser you use:<br>

#### **Per base sequence quality**

|Forward read|Reverse read|
|-|-|
|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20base%20sequence%20quality%20R1.png" width="100%">|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20base%20sequence%20quality%20R2.png" width="100%">|

#### **Per sequence quality scores**

|Forward read|Reverse read|
|-|-|
|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20sequence%20quality%20scores%20R1.png" width="100%">|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20sequence%20quality%20scores%20R2.png" width="100%">|

#### **Per base sequence content**

|Forward read|Reverse read|
|-|-|
|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20base%20sequence%20content%20R1.png" width="100%">|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20base%20sequence%20content%20R2.png" width="100%">|

#### **Per sequence GC content**

|Forward read|Reverse read|
|-|-|
|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20sequence%20GC%20content%20R1.png" width="100%">|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Per%20sequence%20GC%20content%20R2.png" width="100%">|

#### **Adapter Content**

|Forward read|Reverse read|
|-|-|
|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Adapter%20Content%20R1.png" width="100%">|<img src="/images/ngs-handbook/02_Genomic_Variation_Analysis/Adapter%20Content%20R2.png" width="100%">|

The main problem is `Per base sequence quality`. So we will solve it using `trimmomatic`.<br>
Details on `trimmomatic` and QC in the [Quality Control chapter](../01_Quality_Control).

----------------------------------------------

### **Step 3: Trim low-quality bases**

**_Input_**

```bash
mkdir data/trimmed/
```

**_Input_**

```bash
trimmomatic PE -threads 2 data/SRR17909485_1.fastq.gz data/SRR17909485_2.fastq.gz \
    data/trimmed/SRR17909485_R1.trim.paired.fastq.gz data/trimmed/SRR17909485_R1.trim.unpaired.fastq.gz \
    data/trimmed/SRR17909485_R2.trim.paired.fastq.gz data/trimmed/SRR17909485_R2.trim.unpaired.fastq.gz \
    LEADING:22 TRAILING:22 SLIDINGWINDOW:6:22 MINLEN:32
```

**_Output_**

```
TrimmomaticPE: Started with arguments:
 -threads 2 data/SRR17909485_1.fastq.gz data/SRR17909485_2.fastq.gz data/trimmed/SRR17909485_R1.trim.paired.fastq.gz data/trimmed/SRR17909485_R1.trim.unpaired.fastq.gz data/trimmed/SRR17909485_R2.trim.paired.fastq.gz data/trimmed/SRR17909485_R2.trim.unpaired.fastq.gz LEADING:22 TRAILING:22 SLIDINGWINDOW:6:22 MINLEN:32
Quality encoding detected as phred33
Input Read Pairs: 1251776 Both Surviving: 1205173 (96.28%) Forward Only Surviving: 34920 (2.79%) Reverse Only Surviving: 2954 (0.24%) Dropped: 8729 (0.70%)
TrimmomaticPE: Completed successfully
```

- `LEADING:22`: Removes low-quality bases from the start of the reads if they have a quality score lower than 22.
- `TRAILING:22`: Removes low-quality bases from the end of the reads if they have a quality score lower than 22.
- `SLIDINGWINDOW:6:22`: Performs a sliding window trimming. It scans the read with a 6-base wide window, cutting when the average quality per base drops below 22.
- `MINLEN:32`: Discards reads that are shorter than 32 bases after trimming.

Input Read Pairs: 1251776<br>
- Both Surviving: 1205173 (96.28%)<br>
- Forward Only Surviving: 34920 (2.79%)<br>
- Reverse Only Surviving: 2954 (0.24%)<br>
- Dropped: 8729 (0.70%)<br>

----------------------------------------------

### **Step 4: Index reference**

**_Input_**

```bash
bwa index data/reference/EcoliK12MG1655.fa
```

----------------------------------------------

### **Step 5: Align to reference**

**_Input_**

```bash
mkdir data/bam/
```

**_Input_**

```bash
bwa mem -t 2 -R '@RG\tID:1' \
    data/reference/EcoliK12MG1655.fa \
    data/trimmed/SRR17909485_R1.trim.paired.fastq.gz \
    data/trimmed/SRR17909485_R2.trim.paired.fastq.gz \
    | samtools view -b > data/bam/EcoliK12MG1655.SRR17909485.unsorted.bam
```

----------------------------------------------

### **Step 6: Sort alignment**

**_Input_**

```bash
samtools sort --threads 2 data/bam/EcoliK12MG1655.SRR17909485.unsorted.bam > \
    data/bam/EcoliK12MG1655.SRR17909485.sorted.bam
```

----------------------------------------------

### **Step 7: Make bam index**

**_Input_**

```bash
samtools index data/bam/EcoliK12MG1655.SRR17909485.sorted.bam
```

----------------------------------------------

### **Step 8: Realign indels**

**_Input_**

```bash
abra2 --threads 2 --mad 100 --mbq 24 \
    --ref data/reference/EcoliK12MG1655.fa \
    --in data/bam/EcoliK12MG1655.SRR17909485.sorted.bam \
    --out data/bam/EcoliK12MG1655.SRR17909485.final.bam
```

----------------------------------------------

### **Step 9: Index final bam**

**_Input_**

```bash
samtools index data/bam/EcoliK12MG1655.SRR17909485.final.bam
```

----------------------------------------------

### **Step 10: Call variants**

**_Input_**

```bash
mkdir data/vcf/
```

**_Input_**

```bash
bcftools mpileup -Ou --max-depth 5000 \
    -f data/reference/EcoliK12MG1655.fa \
    data/bam/EcoliK12MG1655.SRR17909485.final.bam \
    | bcftools call -mv --ploidy 1 -Ov \
    -o data/vcf/EcoliK12MG1655.SRR17909485.called.bcftools.vcf
```

`Mpileup` is actually a slightly different representation of alignment. `Mpileup` takes very simple things, it takes the number of reads aligned to each point and for each position in the reference genome it counts statistics - how many letters match the reference, how many mismatch, how many insertions, deletions.<br>

Unfortunately `bcftools` is not multithreaded. If we do re-alignment with `abra2`, `bcftools` will work a bit faster.<br>

`bcftools --ploidy 1` - we tell the program that we have a haploid genome.<br>

`vcf` = variant call format.
This is actually a list of variations in a clear text format: chromosome, position, which letter was replaced by which letter, the quality of how much `bcftools` considers it important, some set of certain statistics (frequency of letter occurrence, etc.).<br>

**_Input_**

```bash
tail data/vcf/EcoliK12MG1655.SRR17909485.called.bcftools.vcf
```

**_Output_**

```
U00096.3	3110421	.	C	A	225.417	.	DP=41;VDB=0.454064;SGB=-0.692914;MQSBZ=1.45774;MQ0F=0;AC=1;AN=1;DP4=0,0,8,17;MQ=59	GT:PL	1:255,0
U00096.3	3560455	.	C	CG	225.417	.	INDEL;IDV=82;IMF=1;DP=82;VDB=0.0327919;SGB=-0.693147;MQSBZ=0;BQBZ=0.885506;MQ0F=0;AC=1;AN=1;DP4=0,0,37,45;MQ=60	GT:PL	1:255,0
U00096.3	4093770	.	C	T	225.417	.	DP=141;VDB=0.961761;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,39,52;MQ=60	GT:PL	1:255,0
U00096.3	4161248	.	G	T	225.417	.	DP=117;VDB=0.400301;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,36,33;MQ=60	GT:PL	1:255,0
U00096.3	4164123	.	C	G	225.417	.	DP=130;VDB=0.322528;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,32,44;MQ=60	GT:PL	1:255,0
U00096.3	4296380	.	AC	ACGC	228.422	.	INDEL;IDV=98;IMF=0.989899;DP=99;VDB=8.47905e-05;SGB=-0.693147;RPBZ=-1.71529;MQBZ=1.83497;MQSBZ=-2.15139;BQBZ=-1.59515;SCBZ=-9.8995;MQ0F=0;AC=1;AN=1;DP4=0,1,47,51;MQ=54	GT:PL	1:255,0
U00096.3	4474834	.	A	G	225.417	.	DP=81;VDB=0.345282;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,21,26;MQ=60	GT:PL	1:255,0
U00096.3	4585480	.	G	A	225.417	.	DP=70;VDB=0.996826;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,20,26;MQ=60	GT:PL	1:255,0
U00096.3	4602509	.	C	T	225.417	.	DP=136;VDB=0.843259;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,42,41;MQ=60	GT:PL	1:255,0
U00096.3	4616669	.	G	T	225.417	.	DP=88;VDB=0.58582;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,25,31;MQ=60	GT:PL	1:255,0
```

----------------------------------------------

### **Step 11: Make snpEff database**

`snpEff` has a database of already well-sequenced, annotated genomes of many known organisms. `snpEff` loads the annotations and while annotating the `vcf` file, it checks which gene each particular variation falls into which position.<br>
By _E. coli_ snpEff has over 3500 databases, which is a lot.<br>
The best way to handle this is to make your own database.<br>

**_Input_**

```bash
mkdir -p data/EcoliK12MG1655/
```

**_Input_**

```bash
cp data/reference/EcoliK12MG1655.gb data/EcoliK12MG1655/genes.gbk
```

**_Input_**

```bash
echo "EcoliK12MG1655.genome : EcoliK12MG1655\nEcoliK12MG1655.chromosomes : EcoliK12MG1655.gb\nEcoliK12MG1655.codonTable : Standard" \
    > snpEff.config
```

**_Input_**

```bash
cat snpEff.config
```

**_Output_**

```
EcoliK12MG1655.genome : EcoliK12MG1655
EcoliK12MG1655.chromosomes : EcoliK12MG1655.gb
EcoliK12MG1655.codonTable : Standard
```

1. Line 1 - under what name snpEff will know this reference genome (`EcoliK12MG1655`)
2. Line 2 - the chromosomes from this genome (all contigs with all annotations) will be in the file `EcoliK12MG1655.gb`
3. Line 3 - standard genetic code table

**_Input_**

```bash
snpEff build -c snpEff.config -genbank EcoliK12MG1655
```

**_Output_**

```
00:00:00 Codon table 'Standard' for genome 'EcoliK12MG1655'
	Protein check:	EcoliK12MG1655	OK: 4298	Not found: 17	Errors: 0	Error percentage: 0.0%
```

----------------------------------------------

### **Step 12: Annotate variants**

**_Input_**

```bash
snpEff ann -v EcoliK12MG1655  \
    data/vcf/EcoliK12MG1655.SRR17909485.called.bcftools.vcf > \
    data/vcf/EcoliK12MG1655.SRR17909485.annotated.vcf
```

**_Input_**

```bash
mv snpEff_genes.txt EcoliK12MG1655.SRR17909485.snpEff_genes.txt
```

**_Input_**

```bash
mv snpEff_summary.html EcoliK12MG1655.SRR17909485.snpEff_summary.html
```

**_Input_**

```bash
tail data/vcf/EcoliK12MG1655.SRR17909485.annotated.vcf
```

**_Output_**

```
U00096.3	3110421	.	C	A	225.417	.	DP=41;VDB=0.454064;SGB=-0.692914;MQSBZ=1.45774;MQ0F=0;AC=1;AN=1;DP4=0,0,8,17;MQ=59;ANN=A|upstream_gene_variant|MODIFIER|speC|b2965|transcript|b2965|protein_coding||c.-1266G>T|||||1266|,A|upstream_gene_variant|MODIFIER|yqgH|b4785|transcript|b4785|protein_coding||c.-1021G>T|||||1021|,A|upstream_gene_variant|MODIFIER|yqhJ|b4786|transcript|b4786|protein_coding||c.-845C>A|||||845|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|mltC|b2963|transcript|b2963|protein_coding||c.*4909C>A|||||4909|,A|downstream_gene_variant|MODIFIER|nupG|b2964|transcript|b2964|protein_coding||c.*3451C>A|||||3451|,A|downstream_gene_variant|MODIFIER|yqgA|b2966|transcript|b2966|protein_coding||c.*161C>A|||||161|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|yghD|b2968|transcript|b2968|protein_coding||c.*169G>T|||||169|,A|downstream_gene_variant|MODIFIER|yghE|b2969|transcript|b2969|protein_coding||c.*707G>T|||||707|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|yghF|b2970|transcript|b2970|protein_coding||c.*1633G>T|||||1633|WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|yghG|b2971|transcript|b2971|protein_coding||c.*2646G>T|||||2646|,A|downstream_gene_variant|MODIFIER|pppA|b2972|transcript|b2972|protein_coding||c.*3122G>T|||||3122|,A|downstream_gene_variant|MODIFIER|yghJ|b4466|transcript|b4466|protein_coding||c.*4129G>T|||||4129|,A|intragenic_variant|MODIFIER|pheV|b2967|gene_variant|b2967|||n.3110421C>A||||||	GT:PL	1:255,0
U00096.3	3560455	.	C	CG	225.417	.	INDEL;IDV=82;IMF=1;DP=82;VDB=0.0327919;SGB=-0.693147;MQSBZ=0;BQBZ=0.885506;MQ0F=0;AC=1;AN=1;DP4=0,0,37,45;MQ=60;ANN=CG|upstream_gene_variant|MODIFIER|rtcA|b4475|transcript|b4475|protein_coding||c.-3607_-3606insC|||||3607|,CG|upstream_gene_variant|MODIFIER|rtcB|b3421|transcript|b3421|protein_coding||c.-2377_-2376insC|||||2377|,CG|upstream_gene_variant|MODIFIER|glpD|b3426|transcript|b3426|protein_coding||c.-1558_-1557insG|||||1557|,CG|downstream_gene_variant|MODIFIER|malT|b3418|transcript|b3418|protein_coding||c.*4665_*4666insG|||||4666|,CG|downstream_gene_variant|MODIFIER|rtcR|b3422|transcript|b3422|protein_coding||c.*589_*590insG|||||590|,CG|downstream_gene_variant|MODIFIER|glpG|b3424|transcript|b3424|protein_coding||c.*166_*167insC|||||166|,CG|downstream_gene_variant|MODIFIER|glpE|b3425|transcript|b3425|protein_coding||c.*1041_*1042insC|||||1041|,CG|downstream_gene_variant|MODIFIER|yzgL|b3427|transcript|b3427|protein_coding||c.*3268_*3269insC|||||3268|,CG|downstream_gene_variant|MODIFIER|glgP|b3428|transcript|b3428|protein_coding||c.*3678_*3679insC|||||3678|,CG|intragenic_variant|MODIFIER|glpR|b3423|gene_variant|b3423|||n.3560455_3560456insG||||||	GT:PL	1:255,0
U00096.3	4093770	.	C	T	225.417	.	DP=141;VDB=0.961761;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,39,52;MQ=60;ANN=T|missense_variant|MODERATE|rhaD|b3902|transcript|b3902|protein_coding|1/1|c.503G>A|p.Gly168Asp|503/825|503/825|168/274||,T|upstream_gene_variant|MODIFIER|frvR|b3897|transcript|b3897|protein_coding||c.-3915G>A|||||3915|,T|upstream_gene_variant|MODIFIER|frvX|b3898|transcript|b3898|protein_coding||c.-2845G>A|||||2845|,T|upstream_gene_variant|MODIFIER|frvB|b3899|transcript|b3899|protein_coding||c.-1404G>A|||||1404|,T|upstream_gene_variant|MODIFIER|frvA|b3900|transcript|b3900|protein_coding||c.-947G>A|||||947|,T|upstream_gene_variant|MODIFIER|rhaM|b3901|transcript|b3901|protein_coding||c.-332G>A|||||332|,T|upstream_gene_variant|MODIFIER|rhaS|b3905|transcript|b3905|protein_coding||c.-3966C>T|||||3966|,T|upstream_gene_variant|MODIFIER|rhaR|b3906|transcript|b3906|protein_coding||c.-4876C>T|||||4876|WARNING_TRANSCRIPT_NO_START_CODON,T|downstream_gene_variant|MODIFIER|rhaA|b3903|transcript|b3903|protein_coding||c.*953G>A|||||953|,T|downstream_gene_variant|MODIFIER|rhaB|b3904|transcript|b3904|protein_coding||c.*2209G>A|||||2209|	GT:PL	1:255,0
U00096.3	4161248	.	G	T	225.417	.	DP=117;VDB=0.400301;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,36,33;MQ=60;ANN=T|missense_variant|MODERATE|fabR|b3963|transcript|b3963|protein_coding|1/1|c.182G>T|p.Gly61Val|182/705|182/705|61/234||,T|upstream_gene_variant|MODIFIER|sthA|b3962|transcript|b3962|protein_coding||c.-458C>A|||||458|,T|upstream_gene_variant|MODIFIER|yijD|b3964|transcript|b3964|protein_coding||c.-523G>T|||||523|,T|upstream_gene_variant|MODIFIER|btuB|b3966|transcript|b3966|protein_coding||c.-2391G>T|||||2391|,T|upstream_gene_variant|MODIFIER|murI|b3967|transcript|b3967|protein_coding||c.-4180G>T|||||4180|,T|downstream_gene_variant|MODIFIER|argB|b3959|transcript|b3959|protein_coding||c.*4459G>T|||||4459|,T|downstream_gene_variant|MODIFIER|argH|b3960|transcript|b3960|protein_coding||c.*3025G>T|||||3025|,T|downstream_gene_variant|MODIFIER|oxyR|b3961|transcript|b3961|protein_coding||c.*1841G>T|||||1841|,T|downstream_gene_variant|MODIFIER|trmA|b3965|transcript|b3965|protein_coding||c.*922C>A|||||922|	GT:PL	1:255,0
U00096.3	4164123	.	C	G	225.417	.	DP=130;VDB=0.322528;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,32,44;MQ=60;ANN=G|missense_variant|MODERATE|btuB|b3966|transcript|b3966|protein_coding|1/1|c.485C>G|p.Ala162Gly|485/1845|485/1845|162/614||,G|upstream_gene_variant|MODIFIER|sthA|b3962|transcript|b3962|protein_coding||c.-3333G>C|||||3333|,G|upstream_gene_variant|MODIFIER|trmA|b3965|transcript|b3965|protein_coding||c.-853G>C|||||853|,G|upstream_gene_variant|MODIFIER|murI|b3967|transcript|b3967|protein_coding||c.-1305C>G|||||1305|,G|downstream_gene_variant|MODIFIER|oxyR|b3961|transcript|b3961|protein_coding||c.*4716C>G|||||4716|,G|downstream_gene_variant|MODIFIER|fabR|b3963|transcript|b3963|protein_coding||c.*2352C>G|||||2352|,G|downstream_gene_variant|MODIFIER|yijD|b3964|transcript|b3964|protein_coding||c.*1993C>G|||||1993|	GT:PL	1:255,0
U00096.3	4296380	.	AC	ACGC	228.422	.	INDEL;IDV=98;IMF=0.989899;DP=99;VDB=8.47905e-05;SGB=-0.693147;RPBZ=-1.71529;MQBZ=1.83497;MQSBZ=-2.15139;BQBZ=-1.59515;SCBZ=-9.8995;MQ0F=0;AC=1;AN=1;DP4=0,1,47,51;MQ=54;ANN=ACGC|downstream_gene_variant|MODIFIER|nrfD|b4073|transcript|b4073|protein_coding||c.*4949_*4950insGC|||||4950|,ACGC|downstream_gene_variant|MODIFIER|nrfE|b4074|transcript|b4074|protein_coding||c.*3211_*3212insGC|||||3212|WARNING_TRANSCRIPT_NO_START_CODON,ACGC|downstream_gene_variant|MODIFIER|nrfF|b4075|transcript|b4075|protein_coding||c.*2835_*2836insGC|||||2836|,ACGC|downstream_gene_variant|MODIFIER|nrfG|b4076|transcript|b4076|protein_coding||c.*2242_*2243insGC|||||2243|,ACGC|downstream_gene_variant|MODIFIER|gltP|b4077|transcript|b4077|protein_coding||c.*587_*588insGC|||||588|,ACGC|downstream_gene_variant|MODIFIER|yjcO|b4078|transcript|b4078|protein_coding||c.*54_*55insGC|||||54|,ACGC|downstream_gene_variant|MODIFIER|fdhF|b4079|transcript|b4079|protein_coding||c.*837_*838insGC|||||837|,ACGC|downstream_gene_variant|MODIFIER|mdtP|b4080|transcript|b4080|protein_coding||c.*3182_*3183insGC|||||3182|,ACGC|downstream_gene_variant|MODIFIER|mdtO|b4081|transcript|b4081|protein_coding||c.*4645_*4646insGC|||||4645|,ACGC|intergenic_region|MODIFIER|gltP-yjcO|b4077-b4078|intergenic_region|b4077-b4078|||n.4296381_4296382insGC||||||	GT:PL	1:255,0
U00096.3	4474834	.	A	G	225.417	.	DP=81;VDB=0.345282;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,21,26;MQ=60;ANN=G|upstream_gene_variant|MODIFIER|ridA|b4243|transcript|b4243|protein_coding||c.-3921T>C|||||3921|,G|upstream_gene_variant|MODIFIER|pyrI|b4244|transcript|b4244|protein_coding||c.-3387T>C|||||3387|,G|upstream_gene_variant|MODIFIER|pyrB|b4245|transcript|b4245|protein_coding||c.-2439T>C|||||2439|,G|upstream_gene_variant|MODIFIER|pyrL|b4246|transcript|b4246|protein_coding||c.-2301T>C|||||2301|,G|upstream_gene_variant|MODIFIER|yjgH|b4248|transcript|b4248|protein_coding||c.-1625T>C|||||1625|,G|upstream_gene_variant|MODIFIER|bdcA|b4249|transcript|b4249|protein_coding||c.-781T>C|||||781|,G|upstream_gene_variant|MODIFIER|tabA|b4252|transcript|b4252|protein_coding||c.-28A>G|||||28|,G|upstream_gene_variant|MODIFIER|yjgL|b4253|transcript|b4253|protein_coding||c.-603A>G|||||603|,G|upstream_gene_variant|MODIFIER|rraB|b4255|transcript|b4255|protein_coding||c.-3639A>G|||||3639|,G|upstream_gene_variant|MODIFIER|yjgN|b4257|transcript|b4257|protein_coding||c.-4896A>G|||||4896|,G|downstream_gene_variant|MODIFIER|mgtA|b4242|transcript|b4242|protein_coding||c.*4513A>G|||||4513|,G|downstream_gene_variant|MODIFIER|bdcR|b4251|transcript|b4251|protein_coding||c.*117A>G|||||117|,G|downstream_gene_variant|MODIFIER|argI|b4254|transcript|b4254|protein_coding||c.*2473T>C|||||2473|,G|downstream_gene_variant|MODIFIER|yjgM|b4256|transcript|b4256|protein_coding||c.*4200T>C|||||4200|,G|intergenic_region|MODIFIER|bdcR-tabA|b4251-b4252|intergenic_region|b4251-b4252|||n.4474834A>G||||||	GT:PL	1:255,0
U00096.3	4585480	.	G	A	225.417	.	DP=70;VDB=0.996826;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,20,26;MQ=60;ANN=A|stop_gained|HIGH|hsdR|b4350|transcript|b4350|protein_coding|1/1|c.1282C>T|p.Gln428*|1282/3513|1282/3513|428/1170||,A|upstream_gene_variant|MODIFIER|hsdS|b4348|transcript|b4348|protein_coding||c.-4018C>T|||||4018|,A|upstream_gene_variant|MODIFIER|hsdM|b4349|transcript|b4349|protein_coding||c.-2432C>T|||||2432|,A|upstream_gene_variant|MODIFIER|mrr|b4351|transcript|b4351|protein_coding||c.-1469G>A|||||1469|,A|downstream_gene_variant|MODIFIER|yjiA|b4352|transcript|b4352|protein_coding||c.*2429C>T|||||2429|,A|downstream_gene_variant|MODIFIER|yjiX|b4353|transcript|b4353|protein_coding||c.*3396C>T|||||3396|,A|downstream_gene_variant|MODIFIER|btsT|b4354|transcript|b4354|protein_coding||c.*3649C>T|||||3649|	GT:PL	1:255,0
U00096.3	4602509	.	C	T	225.417	.	DP=136;VDB=0.843259;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,42,41;MQ=60;ANN=T|stop_gained|HIGH|yjjP|b4364|transcript|b4364|protein_coding|1/1|c.350G>A|p.Trp117*|350/771|350/771|117/256||,T|upstream_gene_variant|MODIFIER|opgB|b4359|transcript|b4359|protein_coding||c.-3068G>A|||||3068|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|yjjA|b4360|transcript|b4360|protein_coding||c.-2320G>A|||||2320|,T|upstream_gene_variant|MODIFIER|dnaC|b4361|transcript|b4361|protein_coding||c.-1534G>A|||||1534|,T|upstream_gene_variant|MODIFIER|dnaT|b4362|transcript|b4362|protein_coding||c.-992G>A|||||992|,T|upstream_gene_variant|MODIFIER|yjjB|b4363|transcript|b4363|protein_coding||c.-412G>A|||||412|,T|upstream_gene_variant|MODIFIER|yjjQ|b4365|transcript|b4365|protein_coding||c.-968C>T|||||968|,T|upstream_gene_variant|MODIFIER|bglJ|b4366|transcript|b4366|protein_coding||c.-1651C>T|||||1651|,T|upstream_gene_variant|MODIFIER|yjjZ|b4567|transcript|b4567|protein_coding||c.-3295C>T|||||3295|,T|downstream_gene_variant|MODIFIER|fhuF|b4367|transcript|b4367|protein_coding||c.*2366G>A|||||2366|,T|downstream_gene_variant|MODIFIER|rsmC|b4371|transcript|b4371|protein_coding||c.*4160G>A|||||4160|	GT:PL	1:255,0
U00096.3	4616669	.	G	T	225.417	.	DP=88;VDB=0.58582;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,25,31;MQ=60;ANN=T|missense_variant|MODERATE|yjjI|b4380|transcript|b4380|protein_coding|1/1|c.397C>A|p.Leu133Ile|397/1551|397/1551|133/516||,T|upstream_gene_variant|MODIFIER|yjjW|b4379|transcript|b4379|protein_coding||c.-1126C>A|||||1126|,T|upstream_gene_variant|MODIFIER|deoC|b4381|transcript|b4381|protein_coding||c.-654G>T|||||654|,T|upstream_gene_variant|MODIFIER|deoA|b4382|transcript|b4382|protein_coding||c.-1560G>T|||||1560|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|deoB|b4383|transcript|b4383|protein_coding||c.-2934G>T|||||2934|,T|upstream_gene_variant|MODIFIER|deoD|b4384|transcript|b4384|protein_coding||c.-4214G>T|||||4214|,T|downstream_gene_variant|MODIFIER|osmY|b4376|transcript|b4376|protein_coding||c.*4668G>T|||||4668|,T|downstream_gene_variant|MODIFIER|ytjA|b4568|transcript|b4568|protein_coding||c.*4380G>T|||||4380|,T|downstream_gene_variant|MODIFIER|yjjU|b4377|transcript|b4377|protein_coding||c.*3185G>T|||||3185|WARNING_TRANSCRIPT_NO_START_CODON,T|downstream_gene_variant|MODIFIER|yjjV|b4378|transcript|b4378|protein_coding||c.*2409G>T|||||2409|WARNING_TRANSCRIPT_NO_START_CODON	GT:PL	1:255,0
```

The `vcf` file became significantly larger because the `ANN` (annotation) line appeared after all the technical information<br>
Example:<br>

```
ANN=T|missense_variant|MODERATE|yjjI|b4380|transcript|b4380|protein_coding|1/1|c.397C>A|p.Leu133Ile.
```

- There's been a substitution of a letter for the `T`.
- This is a `missence variant` mutation.
- Its effect is `MODERATE`.
- `yjjI` - gene name
- `b4380` - number of the transcript in which the substitution occurred.
- `protein_coding` - this is the protein coding site
- `c.397C>A` - replacement in nucleotide sequence
- `p.Leu133Ile` - subsequent substitution in protein sequence (133 leucine is replaced by isoleucine)

----------------------------------------------

### **Step 13: Filter high and moderate effect mutations**
Let's filter out high and moderate effect mutations:<br>
- nonsense<br>
- missense<br>
- frameshifts<br>

**_Input_**

```bash
SnpSift filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')" data/vcf/EcoliK12MG1655.SRR17909485.annotated.vcf > data/vcf/EcoliK12MG1655.SRR17909485.higheffect.vcf
```

**_Input_**

```bash
tail data/vcf/EcoliK12MG1655.SRR17909485.higheffect.vcf
```

**_Output_**

```
U00096.3	705013	.	T	C	225.417	.	DP=103;VDB=0.0387749;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,30,36;MQ=60;ANN=C|missense_variant|MODERATE|nagE|b0679|transcript|b0679|protein_coding|1/1|c.1070T>C|p.Leu357Ser|1070/1947|1070/1947|357/648||,C|upstream_gene_variant|MODIFIER|umpH|b0675|transcript|b0675|protein_coding||c.-4687A>G|||||4687|,C|upstream_gene_variant|MODIFIER|nagC|b0676|transcript|b0676|protein_coding||c.-3419A>G|||||3419|,C|upstream_gene_variant|MODIFIER|nagA|b0677|transcript|b0677|protein_coding||c.-2262A>G|||||2262|,C|upstream_gene_variant|MODIFIER|nagB|b0678|transcript|b0678|protein_coding||c.-1402A>G|||||1402|,C|upstream_gene_variant|MODIFIER|glnS|b0680|transcript|b0680|protein_coding||c.-1080T>C|||||1080|,C|upstream_gene_variant|MODIFIER|chiP|b0681|transcript|b0681|protein_coding||c.-3321T>C|||||3321|,C|upstream_gene_variant|MODIFIER|chiQ|b0682|transcript|b0682|protein_coding||c.-4777T>C|||||4777|	GT:PL	1:255,0
U00096.3	1337394	.	A	G	225.417	.	DP=65;VDB=0.915279;SGB=-0.693146;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,24,19;MQ=60;ANN=G|missense_variant|MODERATE|acnA|b1276|transcript|b1276|protein_coding|1/1|c.1564A>G|p.Ser522Gly|1564/2676|1564/2676|522/891||,G|upstream_gene_variant|MODIFIER|pgpB|b1278|transcript|b1278|protein_coding||c.-1936A>G|||||1936|,G|upstream_gene_variant|MODIFIER|lapA|b1279|transcript|b1279|protein_coding||c.-2849A>G|||||2849|WARNING_TRANSCRIPT_NO_START_CODON,G|upstream_gene_variant|MODIFIER|lapB|b1280|transcript|b1280|protein_coding||c.-3164A>G|||||3164|,G|upstream_gene_variant|MODIFIER|pyrF|b1281|transcript|b1281|protein_coding||c.-4527A>G|||||4527|,G|downstream_gene_variant|MODIFIER|topA|b1274|transcript|b1274|protein_coding||c.*3749A>G|||||3749|,G|downstream_gene_variant|MODIFIER|cysB|b1275|transcript|b1275|protein_coding||c.*2565A>G|||||2565|,G|downstream_gene_variant|MODIFIER|ymiA|b4522|transcript|b4522|protein_coding||c.*2106A>G|||||2106|,G|downstream_gene_variant|MODIFIER|yciX|b4523|transcript|b4523|protein_coding||c.*1936A>G|||||1936|,G|downstream_gene_variant|MODIFIER|ymiC|b4741|transcript|b4741|protein_coding||c.*1727A>G|||||1727|,G|downstream_gene_variant|MODIFIER|ribA|b1277|transcript|b1277|protein_coding||c.*1176T>C|||||1176|	GT:PL	1:255,0
U00096.3	1514951	.	G	T	225.417	.	DP=68;VDB=0.981287;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,18,30;MQ=60;ANN=T|missense_variant|MODERATE|ydcV|b1443|transcript|b1443|protein_coding|1/1|c.190G>T|p.Ala64Ser|190/795|190/795|64/264||,T|upstream_gene_variant|MODIFIER|patD|b1444|transcript|b1444|protein_coding||c.-627G>T|||||627|,T|upstream_gene_variant|MODIFIER|ortT|b1445|transcript|b1445|protein_coding||c.-2438G>T|||||2438|,T|upstream_gene_variant|MODIFIER|ydcY|b1446|transcript|b1446|protein_coding||c.-2697G>T|||||2697|,T|upstream_gene_variant|MODIFIER|curA|b1449|transcript|b1449|protein_coding||c.-4076G>T|||||4076|,T|downstream_gene_variant|MODIFIER|ydcR|b1439|transcript|b1439|protein_coding||c.*3542G>T|||||3542|,T|downstream_gene_variant|MODIFIER|ydcS|b1440|transcript|b1440|protein_coding||c.*2152G>T|||||2152|,T|downstream_gene_variant|MODIFIER|ydcT|b1441|transcript|b1441|protein_coding||c.*1121G>T|||||1121|,T|downstream_gene_variant|MODIFIER|ydcU|b1442|transcript|b1442|protein_coding||c.*179G>T|||||179|,T|downstream_gene_variant|MODIFIER|yncL|b4598|transcript|b4598|protein_coding||c.*2148C>A|||||2148|,T|downstream_gene_variant|MODIFIER|ydcZ|b1447|transcript|b1447|protein_coding||c.*2931C>A|||||2931|,T|downstream_gene_variant|MODIFIER|mnaT|b1448|transcript|b1448|protein_coding||c.*3377C>A|||||3377|	GT:PL	1:255,0
U00096.3	2405959	.	G	T	225.417	.	DP=116;VDB=0.00680486;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,42,33;MQ=60;ANN=T|missense_variant|MODERATE|lrhA|b2289|transcript|b2289|protein_coding|1/1|c.683C>A|p.Ala228Glu|683/939|683/939|228/312||,T|upstream_gene_variant|MODIFIER|nuoF|b2284|transcript|b2284|protein_coding||c.-4404C>A|||||4404|,T|upstream_gene_variant|MODIFIER|nuoE|b2285|transcript|b2285|protein_coding||c.-3907C>A|||||3907|,T|upstream_gene_variant|MODIFIER|nuoC|b2286|transcript|b2286|protein_coding||c.-2114C>A|||||2114|,T|upstream_gene_variant|MODIFIER|nuoB|b2287|transcript|b2287|protein_coding||c.-1346C>A|||||1346|,T|upstream_gene_variant|MODIFIER|nuoA|b2288|transcript|b2288|protein_coding||c.-887C>A|||||887|,T|upstream_gene_variant|MODIFIER|alaA|b2290|transcript|b2290|protein_coding||c.-1602G>T|||||1602|,T|upstream_gene_variant|MODIFIER|yfbR|b2291|transcript|b2291|protein_coding||c.-2903G>T|||||2903|,T|downstream_gene_variant|MODIFIER|yfbS|b2292|transcript|b2292|protein_coding||c.*3561C>A|||||3561|WARNING_TRANSCRIPT_NO_START_CODON	GT:PL	1:255,0
U00096.3	4093770	.	C	T	225.417	.	DP=141;VDB=0.961761;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,39,52;MQ=60;ANN=T|missense_variant|MODERATE|rhaD|b3902|transcript|b3902|protein_coding|1/1|c.503G>A|p.Gly168Asp|503/825|503/825|168/274||,T|upstream_gene_variant|MODIFIER|frvR|b3897|transcript|b3897|protein_coding||c.-3915G>A|||||3915|,T|upstream_gene_variant|MODIFIER|frvX|b3898|transcript|b3898|protein_coding||c.-2845G>A|||||2845|,T|upstream_gene_variant|MODIFIER|frvB|b3899|transcript|b3899|protein_coding||c.-1404G>A|||||1404|,T|upstream_gene_variant|MODIFIER|frvA|b3900|transcript|b3900|protein_coding||c.-947G>A|||||947|,T|upstream_gene_variant|MODIFIER|rhaM|b3901|transcript|b3901|protein_coding||c.-332G>A|||||332|,T|upstream_gene_variant|MODIFIER|rhaS|b3905|transcript|b3905|protein_coding||c.-3966C>T|||||3966|,T|upstream_gene_variant|MODIFIER|rhaR|b3906|transcript|b3906|protein_coding||c.-4876C>T|||||4876|WARNING_TRANSCRIPT_NO_START_CODON,T|downstream_gene_variant|MODIFIER|rhaA|b3903|transcript|b3903|protein_coding||c.*953G>A|||||953|,T|downstream_gene_variant|MODIFIER|rhaB|b3904|transcript|b3904|protein_coding||c.*2209G>A|||||2209|	GT:PL	1:255,0
U00096.3	4161248	.	G	T	225.417	.	DP=117;VDB=0.400301;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,36,33;MQ=60;ANN=T|missense_variant|MODERATE|fabR|b3963|transcript|b3963|protein_coding|1/1|c.182G>T|p.Gly61Val|182/705|182/705|61/234||,T|upstream_gene_variant|MODIFIER|sthA|b3962|transcript|b3962|protein_coding||c.-458C>A|||||458|,T|upstream_gene_variant|MODIFIER|yijD|b3964|transcript|b3964|protein_coding||c.-523G>T|||||523|,T|upstream_gene_variant|MODIFIER|btuB|b3966|transcript|b3966|protein_coding||c.-2391G>T|||||2391|,T|upstream_gene_variant|MODIFIER|murI|b3967|transcript|b3967|protein_coding||c.-4180G>T|||||4180|,T|downstream_gene_variant|MODIFIER|argB|b3959|transcript|b3959|protein_coding||c.*4459G>T|||||4459|,T|downstream_gene_variant|MODIFIER|argH|b3960|transcript|b3960|protein_coding||c.*3025G>T|||||3025|,T|downstream_gene_variant|MODIFIER|oxyR|b3961|transcript|b3961|protein_coding||c.*1841G>T|||||1841|,T|downstream_gene_variant|MODIFIER|trmA|b3965|transcript|b3965|protein_coding||c.*922C>A|||||922|	GT:PL	1:255,0
U00096.3	4164123	.	C	G	225.417	.	DP=130;VDB=0.322528;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,32,44;MQ=60;ANN=G|missense_variant|MODERATE|btuB|b3966|transcript|b3966|protein_coding|1/1|c.485C>G|p.Ala162Gly|485/1845|485/1845|162/614||,G|upstream_gene_variant|MODIFIER|sthA|b3962|transcript|b3962|protein_coding||c.-3333G>C|||||3333|,G|upstream_gene_variant|MODIFIER|trmA|b3965|transcript|b3965|protein_coding||c.-853G>C|||||853|,G|upstream_gene_variant|MODIFIER|murI|b3967|transcript|b3967|protein_coding||c.-1305C>G|||||1305|,G|downstream_gene_variant|MODIFIER|oxyR|b3961|transcript|b3961|protein_coding||c.*4716C>G|||||4716|,G|downstream_gene_variant|MODIFIER|fabR|b3963|transcript|b3963|protein_coding||c.*2352C>G|||||2352|,G|downstream_gene_variant|MODIFIER|yijD|b3964|transcript|b3964|protein_coding||c.*1993C>G|||||1993|	GT:PL	1:255,0
U00096.3	4585480	.	G	A	225.417	.	DP=70;VDB=0.996826;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,20,26;MQ=60;ANN=A|stop_gained|HIGH|hsdR|b4350|transcript|b4350|protein_coding|1/1|c.1282C>T|p.Gln428*|1282/3513|1282/3513|428/1170||,A|upstream_gene_variant|MODIFIER|hsdS|b4348|transcript|b4348|protein_coding||c.-4018C>T|||||4018|,A|upstream_gene_variant|MODIFIER|hsdM|b4349|transcript|b4349|protein_coding||c.-2432C>T|||||2432|,A|upstream_gene_variant|MODIFIER|mrr|b4351|transcript|b4351|protein_coding||c.-1469G>A|||||1469|,A|downstream_gene_variant|MODIFIER|yjiA|b4352|transcript|b4352|protein_coding||c.*2429C>T|||||2429|,A|downstream_gene_variant|MODIFIER|yjiX|b4353|transcript|b4353|protein_coding||c.*3396C>T|||||3396|,A|downstream_gene_variant|MODIFIER|btsT|b4354|transcript|b4354|protein_coding||c.*3649C>T|||||3649|	GT:PL	1:255,0
U00096.3	4602509	.	C	T	225.417	.	DP=136;VDB=0.843259;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,42,41;MQ=60;ANN=T|stop_gained|HIGH|yjjP|b4364|transcript|b4364|protein_coding|1/1|c.350G>A|p.Trp117*|350/771|350/771|117/256||,T|upstream_gene_variant|MODIFIER|opgB|b4359|transcript|b4359|protein_coding||c.-3068G>A|||||3068|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|yjjA|b4360|transcript|b4360|protein_coding||c.-2320G>A|||||2320|,T|upstream_gene_variant|MODIFIER|dnaC|b4361|transcript|b4361|protein_coding||c.-1534G>A|||||1534|,T|upstream_gene_variant|MODIFIER|dnaT|b4362|transcript|b4362|protein_coding||c.-992G>A|||||992|,T|upstream_gene_variant|MODIFIER|yjjB|b4363|transcript|b4363|protein_coding||c.-412G>A|||||412|,T|upstream_gene_variant|MODIFIER|yjjQ|b4365|transcript|b4365|protein_coding||c.-968C>T|||||968|,T|upstream_gene_variant|MODIFIER|bglJ|b4366|transcript|b4366|protein_coding||c.-1651C>T|||||1651|,T|upstream_gene_variant|MODIFIER|yjjZ|b4567|transcript|b4567|protein_coding||c.-3295C>T|||||3295|,T|downstream_gene_variant|MODIFIER|fhuF|b4367|transcript|b4367|protein_coding||c.*2366G>A|||||2366|,T|downstream_gene_variant|MODIFIER|rsmC|b4371|transcript|b4371|protein_coding||c.*4160G>A|||||4160|	GT:PL	1:255,0
U00096.3	4616669	.	G	T	225.417	.	DP=88;VDB=0.58582;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=1;AN=1;DP4=0,0,25,31;MQ=60;ANN=T|missense_variant|MODERATE|yjjI|b4380|transcript|b4380|protein_coding|1/1|c.397C>A|p.Leu133Ile|397/1551|397/1551|133/516||,T|upstream_gene_variant|MODIFIER|yjjW|b4379|transcript|b4379|protein_coding||c.-1126C>A|||||1126|,T|upstream_gene_variant|MODIFIER|deoC|b4381|transcript|b4381|protein_coding||c.-654G>T|||||654|,T|upstream_gene_variant|MODIFIER|deoA|b4382|transcript|b4382|protein_coding||c.-1560G>T|||||1560|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|deoB|b4383|transcript|b4383|protein_coding||c.-2934G>T|||||2934|,T|upstream_gene_variant|MODIFIER|deoD|b4384|transcript|b4384|protein_coding||c.-4214G>T|||||4214|,T|downstream_gene_variant|MODIFIER|osmY|b4376|transcript|b4376|protein_coding||c.*4668G>T|||||4668|,T|downstream_gene_variant|MODIFIER|ytjA|b4568|transcript|b4568|protein_coding||c.*4380G>T|||||4380|,T|downstream_gene_variant|MODIFIER|yjjU|b4377|transcript|b4377|protein_coding||c.*3185G>T|||||3185|WARNING_TRANSCRIPT_NO_START_CODON,T|downstream_gene_variant|MODIFIER|yjjV|b4378|transcript|b4378|protein_coding||c.*2409G>T|||||2409|WARNING_TRANSCRIPT_NO_START_CODON	GT:PL	1:255,0
```