# GLORI_mod
**A modified GLORI data processing pipeline, based on the original one from [liucongcas](https://github.com/Liucongcas/GLORI-tools/)**

This pipeline performed several modifications below based on the original codes:

1. use awk to do AtoG change and reversion, which is faster than python
2. Add comments, reformat log file, keep Important tmp files
3. fix bugs: 
    1. When coverge is < 15 for specifi segment (eg. chrM), pandas throw an error when output *.totalCR.txt, add an if-else to advoid.

**Below is the original readme:**

## Table of content
* Background
* Pre-processing the raw sequencing data
* Installation and Requirement
* Example and Usage
* Maintainers and Contributing
* License

## Background
### GLORI
We developed an absolute m6A quantification method (“GLORI”) that is conceptually similar bisulfite sequencing-based quantification of DNA 5-methylcytosine.
GLORI relies on glyoxal and nitrite-mediated deamination of unmethylated adenosines while keeping m6A intact, thereby achieving specific and efficient m6A detection.

## Pre-processing the raw sequencing data
### Annotations
* GLORI-tools exclusively accepts single-end sequencing reads. Prior to use, it is critical to ensure that the input reads are A-to-G converted reads. Therefore, it is important to carefully verify the orientation of the reads before processing them.
* Before inputting data into GLORI-tools, data cleaning is necessary. This involves removing sequencing adapters, low-quality bases, PCR duplicates based on unique molecular identifiers (UMIs), and finally, removing the UMIs themselves.
### Installation
* trim_galore
* seqkit
* FASTX-Toolkit (including fastx_trimmer)
### Example code
* ```trim_galore -q 20 --stringency 1 -e 0.3 --length {length(5’UMI) +25nt} --path_to_cutadapt {cutadapter} --dont_gzip -o {output_dir1} {inputfile}```
* ```seqkit rmdup -j {Thread} -s -D {output_dupname} {filtered_file} > {filter_file2}```
* ```fastx_trimmer -Q 33 -f {length(5′UMI) +1nt} -i {filter_file2} -o {filter_file3}```

## GLORI-tools

GLORI-tools is a bioinformatics pipeline tailored for the analysis of high-throughput sequencing data generated by GLORI.
![GLORI pipeline](./GLORI-pipeline.jpg " GLORI pipeline ")

## Installation and Requirement
### Installation
GLORI-tools is written in Python3 and is executed from the command line. To install GLORI-tools simply download the code and extract the files into a GLORI-tools installation folder.

### GLORI-tools needs the following tools to be installed and ideally available in the PATH environment:
* STAR ≥ v2.7.5c
* bowtie (bowtie1) version ≥ v1.3.0
* samtools ≥ v1.10
* python ≥ v3.8.3

### GLORI-tools needs the following python package to be installed:
pysam,pandas,argparse,time,collections,os,sys,re,subprocess,multiprocessing,copy,numpy,scipy (v1.10.0),math,sqlite3,Bio,statsmodels,itertools,heapq,glob,signal

## Example and Usage:

### 1. Generate annotation files (required)
#### 1.1 download files for annotation (required, using hg38 as example): 
* ``` wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt ```
* ``` wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz ```

#### 1.2 download reference genome and transcriptome
* ``` wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ```

* ```wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz```

#### 1.3 Unify chromosome naming in GTF file and genome file:
 
* ```python ./get_anno/change_UCSCgtf.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf -j GCF_000001405.39_GRCh38.p13_assembly_report.txt -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens ```

### 2. get reference for reads alignment (required)

#### 2.1 build genome index using STAR

* ``` python ./pipelines/build_genome_index.py -f $genome_fastafile -p 20 -pre hg38 ```

you will get:
* $ hg38.rvsCom.fa
* $ hg38.AG_conversion.fa
* the corresponding index from STAR

#### 2.2 build transcriptome index using bowtie

2.2.1 get the longest transcript for genes (required)

Get the required annotation table files:

* ```python ./get_anno/gtf2anno.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl``` 

* ```awk '$3!~/_/&&$3!="na"' GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl | sed '/unknown_transcript/d'  > GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2```

Get the longest transcript:

* ``` python ./get_anno/selected_longest_transcrpts_fa.py -anno GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2 -fafile GCF_000001405.39_GRCh38.p13_rna.fa --outname_prx GCF_000001405.39_GRCh38.p13_rna2.fa```

2.2.2 build reference with bowtie

* ```python ./pipelines/build_transcriptome_index.py -f $ GCF_000001405.39_GRCh38.p13_rna2.fa -pre GCF_000001405.39_GRCh38.p13_rna2.fa```

you will get:
* $ GCF_000001405.39_GRCh38.p13_rna2.fa.AG_conversion.fa
* the corresponding index from bowtie

### 3. get_base annotation (optional)

#### 3.1 get annotation at single-base resolution

* ```python ./get_anno/anno_to_base.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2 -threads 20 -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno```

#### 3.2 get required annotation file for further removal of duplicated loci

* ``` python ./get_anno/gtf2genelist.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens -f GCF_000001405.39_GRCh38.p13_rna.fa -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist > output2```

* ```awk '$6!~/_/&&$6!="na"' GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist > GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2```

#### 3.3 Removal of duplicated loci in the annotation file

* ```python ./get_anno/anno_to_base_remove_redundance_v1.0.py -i GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno -o GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base -g GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2```

Finally, you will get annotation files: 
* GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base

### 4. alignment and call sites (required)

GLORI-tools takes cleaned reads as input and finally reports files for the conversion rate (A-to-G) of GLORI for each gene and m6A sites at single-base resolution with corresponding A rate representative for modification level. 

#### 4.1 Example shell scripts

| Used files |
| :--- |
| Thread=1 |
| genomdir=your_dir |
| genome=${genomdir}/hg38.AG_conversion.fa |
| genome2=${genomdir}/hg38.fa |
| rvsgenome=${genomdir}/hg38.revCom.fa |
| TfGenome=${genomdir}/GCF_000001405.39_GRCh38.p13_rna2.fa.AG_conversion.fa |
| annodir=your_dir |
| baseanno=${annodir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base |
| anno=${annodir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2 |
| outputdir=your_dir |
| tooldir=/tool_dir/GLORI -tools |
| prx=your_prefix |
| file=your_cleaned_reads | 

#### 4.2 Call m6A sites annotated the with genes

* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b $baseanno -pre ${prx} -o $outputdir --combine --rvs_fac ```

#### 4.3 Call m6A sites without annotated genes.

* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -pre ${prx} -o $outputdir --combine --rvs_fac ```

* In this situation, the background for each m6A sites is the overall conversion rate.

* The site list obtained by the above two methods is basically the similar, and there may be a few differential sites in the list.

#### 4.4 mapping with samples without GLORI treatment

* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -pre ${prx} -o $outputdir --combine --untreated ```

#### 4.5 call all the annotated A cites

1) with annotation:
* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b $baseanno -pre ${prx} -o $outputdir --combine --rvs_fac -c 1 -C 0 -r 0 -p 1.1 -adp 1.1 -s 0 ```
  
2) or without annotation
* ``` python ${tooldir}/run_GLORI.py -i $tooldir -q ${file} -T $Thread -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -pre ${prx} -o $outputdir --combine --rvs_fac -c 1 -C 0 -r 0 -p 1.1 -adp 1.1 -s 0 ```

### 5 Resultes:
#### 5.1 Output files of GLORI

| Output files | Interpretation |
| :---: | :---: |
| ${your_prefix}_merged.sorted.bam | The overall mapping of reads in .bam format |
|  ${your_prefix}_referbase.mpi |  The text pileup output from .bam files |
|  ${your_prefix}.totalCR.txt | The text file containing the median value of the overall A-to-G conversion rate for each transcriptome and gene |
|  ${your_prefix}.totalm6A.FDR.csv | The final list of m6A sites obtained, with the A rate serving as the m6A level |

#### 5.2 GLORI sites files

| Columns | Interpretation |
| :---: | :---: |
| Chr | chromosome |
| Sites | genomic loci |
| Strand | strand |
| Gene | annotated gene |
| CR | conversion rate for genes |
| AGcov | reads coverage with A and G |
| Acov | reads coverage with A |
| Genecov | mean coverage for the whole gene |
| Ratio | A rate for the sites/ or methylation level for the sites |
| Pvalue | test for A rate based on the background |
| P_adjust | FDR ajusted P value |

#### 5.3 GLORI conversion rate (.totalCR) files
| Columns | Interpretation |
| :---: | :---: |
| SA | chromosomes or genes |
| A-to-G_ratio |  A-to-G conversion rate for each chromosome and gene |

## Maintainers and Contributing
* GLORI-tools is developed and maintained by Cong Liu (liucong-1112@pku.edu.cn).
* The development of GLORI-tools is inseparable from the open source software RNA-m5C (https://github.com/SYSU-zhanglab/RNA-m5C).

## Licences
* Released under MIT license

## Citation
please cite:
Liu C., Sun H., Yi Y., Shen W., Li K., Xiao Y., et al. (2022). Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI. Nat. Biotechnol. (https://www.nature.com/articles/s41587-022-01487-9)




