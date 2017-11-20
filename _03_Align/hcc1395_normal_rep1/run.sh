#!/bin/bash
set -e 

STAR --readFilesIn _02_QC/hcc1395_normal_rep1/cleanData/hcc1395_normal_rep1_r1_trimmed.fq.gz _02_QC/hcc1395_normal_rep1/cleanData/hcc1395_normal_rep1_r2_trimmed.fq.gz --outFileNamePrefix _03_Align/hcc1395_normal_rep1/STAR_ --runThreadN 2 --genomeDir genome/referenceFile/ --sjdbGTFfile genome/chr22_with_ERCC92.gtf --genomeLoad NoSharedMemory --outFilterType BySJout --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM MD XS --twopassMode Basic --quantMode GeneCounts --outStd Log --outFilterIntronMotifs RemoveNoncanonical --readFilesCommand zcat --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --outFilterMultimapScoreRange 1 --sjdbScore 2 --sjdbOverhang 100 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.33 --alignSJDBoverhangMin 1 --alignMatesGapMax 1000000 --alignIntronMax 500000 

java -XX:-UsePerfData -Xms8G -Xmx16G -jar $PICARD MarkDuplicates VERBOSITY=WARNING QUIET=false VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2500000 ASSUME_SORTED=false CREATE_INDEX=true REMOVE_DUPLICATES=false INPUT=_03_Align/hcc1395_normal_rep1/accepted_hits.bam OUTPUT=_03_Align/hcc1395_normal_rep1/mk_dup.bam METRICS_FILE=_03_Align/hcc1395_normal_rep1/mkDup_metrics.txt

