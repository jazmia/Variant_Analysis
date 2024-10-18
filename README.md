# Directed Evolution Variant Analysis

This repository contains various scripts for analysing genomic variants of a directed evolution study.

## Installation ## 
Create the conda environment:  `conda create -f environment.yml`.

## Usage ##

**FastQC on sequencing files with FastQC and MultiQC**
```
fastqc sample_seq.fq
multiqc fastq_seq
```

**Trimming of adapter sequences with trimmomatic**
```
trimmomatic PE \
    file1.fq file2.fq \
    file1.trim.fq file1un.trim.fq \
    file2.trim.fq file2un.trim.fq \
    SLIDINGWINDOW:4:20 MINLEN:25 \
    ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
```

**Variant analysis with samtools and bcftools**\
Reference genome: GenBank CP110619.1
```
bwa index reference_seq
bwa mem reference_seq sample_seq.fq > aligned_seq.sam
samtools view -S -b aligned_seq.sam > aligned_seq.bam
samtools sort -o aligned_sorted_seq.bam aligned_seq.bam
samtools flagstat aligned_sorted_seq.bam
bcftools mpileup -O b -o sample_seq.bcf reference_seq akigned_sorted_seq.am
bcftools call --ploidy 1 -m -v -o sample_varaints.vcf sample_seq.bcf
vcfutils.pl varFilter sample_variants.vcf > sample_final_variants.vcf
less -S sample_final_variants.vcf
```

**Generate consensus sequences**\
Done for parental strain used in directed evolution and used in variant analysis of evolved strains\
Also done for samples but not used in variant analysis
```
samtools consensus -f fasta aligned_sorted_seq.bam -o parental.fasta
```

**Annotation with prokka**
```
prokka --outdir annotated_seq --prefix annotated sample_consensus_seq.fasta
```

**Visualisation**
```
samtools index aligned_sorted_seq.bam
download bam bai and vcf files
annotated download fasta and gff files
```

**SNIPPY**\
Do not use consensus genomes from previous step
```
snippy --outdir sample_snps --ref reference_seq --R1 file1.trim.fq --R2 file2.trim.fq
```

# Phenotypic follow-up

Additional R scripts utilised to analyse growth curve data.

## Scripts ##
**Growth_Curve_Analysis_Script**
```
R script to analyse growth curve data including area under curve (AUC) analysis. \
data: 20240913_RIF_fabF.xlsx
      Plate_Layout.xlsx
      20241017_fabF_biolog.xlsx
```
