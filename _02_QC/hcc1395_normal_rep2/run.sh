#!/bin/bash
set -e 

fastqc --threads 2 --extract --quiet --outdir _02_QC/hcc1395_normal_rep2 _01_rawData/hcc1395_normal_rep2_r1.fastq.gz _01_rawData/hcc1395_normal_rep2_r2.fastq.gz

trim_galore --quality 20 --phred33 --fastqc_args "--extract --quiet" --illumina --length 20 --stringency 3 -e 0.1 --output_dir _02_QC/hcc1395_normal_rep2 --paired _01_rawData/hcc1395_normal_rep2_r1.fastq.gz _01_rawData/hcc1395_normal_rep2_r2.fastq.gz
