### <a name="QualChk-1"></a>1. 测序数据质量评估

#### <a name="QualChk-1-1"></a>1.1 测序数据过滤

&nbsp;&nbsp;&nbsp;&nbsp;测序得到的原始下机序列通常会存在测序接头序列以及低质量序列，为了确保分析数据的质量以及结果的可靠性，需要先将这部分序列去除，得到高质量的clean reads，才能用于下一步的分析。  
&nbsp;&nbsp;&nbsp;&nbsp;用trim_galore（version 0.4.3）进行去接头、去低质量处理。Q值（Phred Quality Score）设定为30。对处理的结果进行统计，如下表所示。

**表1.1 数据质量评估表**

|**Sample**	|**Raw_reads**	|**Raw_bases**	|**Clean Reads**	|**Clean Bases**	|**MeanLength**	|**Adaptor(Ratio)**	|
|:---------	|:------------	|:------------	|:--------------	|:--------------	|:-------------	|:-----------------	|
|hcc1395_normal_rep1_r1|331,958|50,125,658 bp|331,958 (100.0%)|48,983,230 bp (97.7%)|147|73,101 (22.0%)|
|hcc1395_normal_rep1_r2|331,958|50,125,658 bp|331,958 (100.0%)|48,195,531 bp (96.1%)|145|71,931 (21.7%)|
|hcc1395_normal_rep2_r1|331,958|50,125,658 bp|331,958 (100.0%)|48,972,517 bp (97.7%)|147|74,846 (22.5%)|
|hcc1395_normal_rep2_r2|331,958|50,125,658 bp|331,958 (100.0%)|48,311,967 bp (96.4%)|145|73,828 (22.2%)|
|hcc1395_tumor_rep1_r1|390,607|58,981,657 bp|390,607 (100.0%)|57,291,836 bp (97.1%)|146|109,003 (27.9%)|
|hcc1395_tumor_rep1_r2|390,607|58,981,657 bp|390,607 (100.0%)|56,506,450 bp (95.8%)|144|105,411 (27.0%)|
|hcc1395_tumor_rep2_r1|390,607|58,981,657 bp|390,607 (100.0%)|57,278,618 bp (97.1%)|146|110,413 (28.3%)|
|hcc1395_tumor_rep2_r2|390,607|58,981,657 bp|390,607 (100.0%)|56,601,236 bp (96.0%)|144|107,563 (27.5%)|

1. *Sample：样品名称*
2. *Raw reads：统计原始序列数据，以四行为一个单位，统计每个文件的测序序列的个数。*
3. *Raw bases：统计原始序列数据，为Raw_reads乘以每条reads的长度。*
4. *Clean reads：计算方法同 Raw Reads，只是统计的文件为过滤后的测序数据。*
5. *Clean bases：Clean reads的个数乘以reads过滤的长度。*
6. *MeanLength：过滤后的reads长度的平均值。*
7. *Adaptor（ratio）：含有接头序列reads占所有reads的比例。*

#### <a name="QualChk-1-2"></a>1.2 测序数据质量评估

##### <a name="QualChk-1-2-1"></a>1.2.1 测序错误分布率检查

&nbsp;&nbsp;&nbsp;&nbsp;测序错误率分布用于检测在测序长度范围内，有无异常的碱基位置存在高错误率。每个碱基测序错误率是测序Phred数值（Phred score, Qphred）通过公式：Qphred = -10log10(e) （e为测序错误率）转化得到，而Phred 数值是在碱基识别（Base Calling）过程中通过一种预测碱基判别发生错误概率模型计算得到的，对应关系如下表：

**表1.2.1 Illumina Casava碱基识别与Phred分值对应关系**

|Phred 分值	|不正确的碱基识别	|碱基正确识别率	|Q-score|
|:--------	|:-------------	|:-----------	|:-----	|
|10	|1/10	|90%	|Q10|
|20	|1/100	|99%	|Q20|
|30	|1/1000	|99.9%	|Q30|
|40	|1/10000	|99.99%	|Q40|

&nbsp;&nbsp;&nbsp;&nbsp;测序错误率与碱基质量有关，受测序仪本身、测序试剂、样品等多个因素共同影响。对于RNA-seq技术，测序错误率通常具有两个特征：1）测序错误率会随着测序序列长度的增加而升高，这是由于测序过程中化学试剂的消耗而导致的。2）前6个碱基的位置也会发生较高的测序错误率，而这个长度正好等于RNA-seq建库过程中反转录所需要的随机引物的长度。因此，推测前6个碱基测序错误率较高是随机引物和RNA模板的不完全结合导致的。

| **Before Quality Trimming**	| **After Quality Trimming**	|
|----------------------------	|----------------------------	|
|*hcc1395_normal_rep1_r1*| |
|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_normal_rep1_r2*| |
|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_normal_rep2_r1*| |
|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_normal_rep2_r2*| |
|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_tumor_rep1_r1*| |
|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_tumor_rep1_r2*| |
|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_tumor_rep2_r1*| |
|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_sequence_quality.png" width="400" height="300">|
|*hcc1395_tumor_rep2_r2*| |
|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_fastqc/Images/per_base_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_base_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_base_quality.png" width="400" height="300">|
|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_sequence_quality.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_sequence_quality.png" width="400" height="300">|

1. *上图左：横坐标为reads碱基坐标，纵坐标为碱基质量值，蓝色线代表平均质量值*
2. *上图右：横坐标为平均序列质量，纵坐标为reads数目*

##### <a name="QualChk-1-2-2"></a>1.2.2 A/T/G/C含量分布检查

&nbsp;&nbsp;&nbsp;&nbsp;正常情况下如果建库均匀，四种碱基的出现频率应该是接近的，而且没有位置差异。GC含量分布检查用于检测有无AT、GC分离现象，而这种现象可能是测序或者建库带来的，并且会影响后续定量分析。而在测序中，由于随机引物扩增偏差等原因，常常会导致在测序得到的每个read前6-7个碱基有较大的波动，这种波动属于正常现象。

| **序列碱基分布图**	| **序列GC含量分布图**	|
|----------------------------	|----------------------------	|
|*hcc1395_normal_rep1_r1*| |
|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r1_val_1_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_normal_rep1_r2*| |
|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_normal_rep1/hcc1395_normal_rep1_r2_val_2_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_normal_rep2_r1*| |
|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r1_val_1_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_normal_rep2_r2*| |
|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_normal_rep2/hcc1395_normal_rep2_r2_val_2_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_tumor_rep1_r1*| |
|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r1_val_1_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_tumor_rep1_r2*| |
|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_tumor_rep1/hcc1395_tumor_rep1_r2_val_2_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_tumor_rep2_r1*| |
|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r1_val_1_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|
|*hcc1395_tumor_rep2_r2*| |
|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_base_sequence_content.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_base_sequence_content.png" width="400" height="300">|<img src="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_sequence_gc_content.png" alt="_02_QC/hcc1395_tumor_rep2/hcc1395_tumor_rep2_r2_val_2_fastqc/Images/per_sequence_gc_content.png" width="400" height="300">|

1. *上图左：横坐标为reads碱基坐标，纵坐标为碱基含量（%）。好的样本中四条线应该平行且接近。当部分位置碱基的比例出现bias时，即四条线在某些位置纷乱交织，往往提示我们有overrepresented sequence的污染*
2. *上图右：横坐标为平均GC含量（%），纵坐标为reads数。红线是实际分布情况，蓝线是理论分布（正态分布）。曲线形状的偏差往往是由于文库的污染或是部分reads构成的子集有偏差（overrepresented reads）。形状接近正态但偏离理论分布的情况提示我们可能有系统偏差*


