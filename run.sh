# read trimming and quality check
qualchk +T template/1.QualChk.md +p "--paired" +t 2 +o "--output_dir _02_QC" +l _01_rawData/filelist.txt

# alignment with STAR
align +i "_02_QC/hcc1395_normal_rep1/cleanData/hcc1395_normal_rep1_r1_trimmed.fq.gz _02_QC/hcc1395_normal_rep1/cleanData/hcc1395_normal_rep1_r2_trimmed.fq.gz" +o _03_Align/STAR +t 2 +g genome/referenceFile/ +G genome/chr22_with_ERCC92.gtf +P STAR ++isMkDup False ++AlnArgs "--outSAMattributes NH HI AS nM NM MD XS --outFilterMatchNminOverLread 0.3 --runThreadN 3"
