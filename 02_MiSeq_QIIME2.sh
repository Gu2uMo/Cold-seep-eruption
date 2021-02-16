## LMO
###----install QIMME2 (2020.8)----
wget https://data.qiime2.org/distro/core/qiime2-2020.8-py36-linux-conda.yml
conda env create -n qiime2-2020.8 --file qiime2-2020.8-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2020.8-py36-linux-conda.yml

###----prepare database for Qiime2  (Silva_138)----
#download the resource from https://docs.qiime2.org/2020.8/data-resources/
source activate qiime2-2020.8
wget https://data.qiime2.org/2020.8/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2020.8/common/silva-138-99-tax.qza
###----train the classifier----
####primer for archaea 516F TGYCAGCCGCCGCGGTAAHACCVGC 855R TCCCCCGCCAATTCCTTTAA
#####extract the corresponding region 
qiime feature-classifier extract-reads \
	--i-sequences ../silva-138-99-seqs.qza \
	--p-f-primer TGYCAGCCGCCGCGGTAAHACCVGC \
	--p-r-primer TCCCCCGCCAATTCCTTTAA \
	--o-reads  00_archaea_516F_855R.qza \
	--p-n-jobs 39 \
	--verbose > log_03_01_arc_train_database.out 
#####train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
	--i-reference-reads 00_archaea_516F_855R.qza \
	--i-reference-taxonomy ../silva-138-99-tax.qza \
	--o-classifier Arc_516F_855R_classifier.qza \
	--verbose > log_03_02_classifier.out 
#####change to the path where the data was stored and link the classifier
#####meta file, saved as 00_sample_info.txt
sample-id       barcode Treat   Group   plot_order      NAME
AOrigin AAAAC   Origin  Origin  A       Origin
A8MPa   AAAAG   8MPa    FC      B       8MPa
A15MPa  AAAAA   15MPa   FC      C       15MPa
A30MPa  AAATG   30MPa   FC      D       30MPa
A8MPaII AAAAT   8MPa    SC      E       8MPaII
A15MPaII        AAAAA   15MPa   SC      F       15MPaII
A30MPaII        AAAAA   30MPa   SC      G       30MPaII
#####load the input file named "00_Cadiz_arc.fa"
qiime tools import --type SampleData[Sequences] \
	--input-path 00_Cadiz_arc.fa \
	--output-path 01_import.qza
qiime vsearch dereplicate-sequences \
    --i-sequences 01_import.qza \
    --o-dereplicated-table 02_vsearch_dep.table.qza \
    --o-dereplicated-sequences 02_vsearch_dep.seq.qza
qiime feature-table summarize \
    --i-table 02_vsearch_dep.table.qza \
    --o-visualization 03_vsearch_table.qzv
#####run the classifier 
qiime feature-classifier classify-sklearn \
    --i-classifier Arc_516F_855R_classifier.qza \
    --p-n-jobs 39 \
    --i-reads 02_vsearch_dep.qza \
    --o-classification 04_vsearch_classification.qza
#####remove the bacteria from the classification file and table file
qiime taxa filter-table --i-table 02_vsearch_dep.table.qza \
    --i-taxonomy 04_vsearch_classification.qza \
    --p-exclude Bacteria,Eukaryota,Unassigned \
    --o-filtered-table 05_taxonomy_filter.table.qza
#####plot the classification, for the abundance matrix
qiime taxa barplot --i-table 05_taxonomy_filter.table.qza \
    --i-taxonomy 04_vsearch_classification.qza \
    --m-metadata-file 00_sample_info.txt \
    --o-visualization 07_tax_bar_plot.qzv
#####calculate the evenness
qiime feature-table summarize \
    --i-table 05_taxonomy_filter.table.qza \
    --o-visualization 05_taxonomy_filter.table.qzv
##### from "05_taxonomy_filter.table.qzv" find the minimum feature frequency
minfreq=32152
qiime diversity alpha-rarefaction \
    --i-table 05_taxonomy_filter.table.qza \
    --p-max-depth ${minfreq} --p-metrics heip_e \
    --p-metrics pielou_e --p-metrics simpson_e \
    --p-metrics chao1 --p-metrics shannon \
    --m-metadata-file 00_sample_info.txt \
    --o-visualization 08_not_tree_alpha_div_${minfreq}


####primer for bacteria 533F TGCCAGCAGCCGCGGTAA 806R GGACTACCAGGGTATCTAATCCTGTT, the work flow is the same as archaea except the primer sequence 
#####extract the corresponding region 
qiime feature-classifier extract-reads \
	--i-sequences ../silva-138-99-seqs.qza \
	--p-f-primer TGCCAGCAGCCGCGGTAA \
	--p-r-primer GGACTACCAGGGTATCTAATCCTGTT \
	--o-reads  00_bacteria_533F_806R.qza \
	--p-n-jobs 39 \
	--verbose > log_03_01_arc_train_database.out 
#####train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
	--i-reference-reads 00_bacteria_533F_806R.qza \
	--i-reference-taxonomy ../silva-138-99-tax.qza \
	--o-classifier Bac_533F_806R_classifier.qza \
	--verbose > log_03_02_classifier.out 

#####meta file, saved as 00_sample_info.txt
sample-id       barcode Treat   Group   plot_order      NAME
BOrigin AAAAC   Origin  Origin  A       Origin
B8MPa   AAAAG   8MPa    FC      B       8MPa
B15MPa  AAAAT   15MPa   FC      C       15MPa
B30MPa  AAAAG   30MPa   FC      D       30MPa
B8MPaII AAAAT   8MPa    SC      E       8MPaII
B15MPaII        AAAAA   15MPa   SC      F       15MPaII
B30MPaII        AAAAG   30MPa   SC      G       30MPaII
#####load the input file named "00_Cadiz_bac.fa"
qiime tools import --type SampleData[Sequences] \
	--input-path 00_Cadiz_bac.fa \
	--output-path 01_import.qza
qiime vsearch dereplicate-sequences \
    --i-sequences 01_import.qza \
    --o-dereplicated-table 02_vsearch_dep.table.qza \
    --o-dereplicated-sequences 02_vsearch_dep.seq.qza
qiime feature-table summarize \
    --i-table 02_vsearch_dep.table.qza \
    --o-visualization 03_vsearch_table.qzv
#####run the classifier 
qiime feature-classifier classify-sklearn \
    --i-classifier Bac_533F_806R_classifier.qza \
    --p-n-jobs 39 \
    --i-reads 02_vsearch_dep.qza \
    --o-classification 04_vsearch_classification.qza
#####remove the bacteria from the classification file and table file
qiime taxa filter-table --i-table 02_vsearch_dep.table.qza \
    --i-taxonomy 04_vsearch_classification.qza \
    --p-exclude Archaea,Eukaryota,Unassigned \
    --o-filtered-table 05_taxonomy_filter.table.qza
#####plot the classification
qiime taxa barplot --i-table 05_taxonomy_filter.table.qza \
    --i-taxonomy 04_vsearch_classification.qza \
    --m-metadata-file 00_sample_info.txt \
    --o-visualization 07_tax_bar_plot.qzv
#####calculate the evenness
qiime feature-table summarize \
    --i-table 05_taxonomy_filter.table.qza \
    --o-visualization 05_taxonomy_filter.table.qzv
##### from "05_taxonomy_filter.table.qzv" find the minimum feature frequency
minfreq=47606
qiime diversity alpha-rarefaction \
    --i-table 05_taxonomy_filter.table.qza \
    --p-max-depth ${minfreq} --p-metrics heip_e \
    --p-metrics pielou_e --p-metrics simpson_e \
    --p-metrics chao1 --p-metrics shannon \
    --m-metadata-file 00_sample_info.txt \
    --o-visualization 08_not_tree_alpha_div_${minfreq}
