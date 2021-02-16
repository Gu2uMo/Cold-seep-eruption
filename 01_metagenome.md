#### 1.raw_data -> clean data

bbduk to remove the low quality reads (Q20, trimmed length <90 bp, contain adapeters)

```bash
sample_list=(Origin L8 L15 L30 L8II H15 H30)
for sample in ${sample_list[@]}
do
		r1=${sample}_00_raw_1.fq.gz
		r2=${sample}_00_raw_2.fq.gz
		c1=${name}_01_clean_bbduk_1.fq.gz
		c2=${name}_01_clean_bbduk_2.fq.gz
		bbduk.sh in=${r1} in2=${r2} out=${c1} out2=${c2} \
        	qtrim=rl trimq=20 mlf=0.6 threads=40 ref=adapters
done
```

tadpole.sh to filter reads containing low-frequency *k*mer

```bash
sample_list=(Origin L8 L15 L30 L8II H15 H30)
for sample in ${sample_list[@]}
do
    f1c=${ele}_01_clean_bbduk_1.fq.gz
    f2c=${ele}_01_clean_bbduk_2.fq.gz
    f1t=${ele}_03_tadpole_1.fq.gz
    f2t=${ele}_03_tadpole_2.fq.gz
    tadpole.sh in=${f1c} in2=${f2c} out=${f1t} out2=${f2t} \
            threads=40 ecc=t tossjunk=t tossdepth=2 \
            tossuncorrectable mbb=2 
done
```

#### 2. assembly

Co-assembly by SPAdes with merged reads. Then scaffolds <500 bp were discarded

```bash
cat 00_reads/*tadpole*_1.fq.gz > 01_clean_1.fq.gz
cat 00_reads/*tadpole*_2.fq.gz > 01_clean_2.fq.gz
spades.py -1 01_clean_1.fq.gz -2 01_clean_2.fq.gz -o 03_spades \
		--meta --only-assembler -t 40 -k 65,75,95,115,127
01_length_cut.py 03_spades/scaffolds.fasta 03_spades.fa 500
```

#### 3. reads mapping and depth calculation

```bash
mkdir 04_map_bbmap
cd 04_map_bbmap
ln -s ../03_spades.fa .
sample_list=(Origin L8 L15 L30 L8II H15 H30)
for ele in ${sample_list[@]}
do
    f1c=../00_reads/${ele}_03_tadpole_1.fq.gz
    f2c=../00_reads/${ele}_03_tadpole_2.fq.gz
    bbmap.sh in=${f1c} in2=${f2c} ref=03_spades.fa nodisk k=13 minid=0.95  \
        keepnames=t covstats=${ele}_04_map_bbmap.depth minaveragequality=5 \
        out=${ele}_04_map_bbmap_temp.sam trimreaddescriptions=t pairlen=350 \
        rescuedist=650
    samtools view -@ 39 -b ${ele}_04_map_bbmap_temp.sam |samtools sort -@ 39 -o ${ele}_04_map_bbmap.bam
    samtools index ${ele}_04_map_bbmap.bam -@ 39
    rm ${ele}_04_map_bbmap_temp.sam 
done 
#script from MetaBat is also used to calculate the depth profile, which is only used for binning
jgi_summarize_bam_contig_depths --outputDepth Cadiz_merge_05_depth_bbmap.metabat \
    ../04_map_bbmap/Cadiz_Origin_04_map_bbmap.bam \
    ../04_map_bbmap/Cadiz_8MPa_04_map_bbmap.bam \
    ../04_map_bbmap/Cadiz_15MPa_04_map_bbmap.bam \
    ../04_map_bbmap/Cadiz_30MPa_04_map_bbmap.bam \
    ../04_map_bbmap/Cadiz_8MPaII_04_map_bbmap.bam \
    ../04_map_bbmap/Cadiz_15MPaII_04_map_bbmap.bam \
    ../04_map_bbmap/Cadiz_30MPaII_04_map_bbmap.bam
```

####  4. Gene annotation abundance and TPM

```bash
refinem call_genes -x fa -c 39 . 01_gene
01_length_cut.py 01_gene/03_spades_genes.faa Cadiz_merge_03_gene_prodigal.faa 33
01_gff_cut.py 01_gene/03_spades_genes.gff Cadiz_merge_03_gene_prodigal.gff 33

#integrate results for KO, GhostKOALA > KofamKOALA > eggnog
##GhostKOALA https://www.kegg.jp/ghostkoala/
##KofamKOALA
~/software/kofamkoala/1.3/kofam_scan-1.3.0/exec_annotation \
        -c ~/software/kofamkoala/1.3/kofam_scan-1.3.0/config-template.yml \
        -o Cadiz_merge_05_ann_kofamkoala.txt \
        -f detail-tsv -E 1e-5 \
        ../03_gene/Cadiz_merge_03_gene_prodigal.faa
##eggnog
emapper.py -o Cadiz_merge_05_ann_eggnog.txt -i Cadiz_merge_03_gene_prodigal.faa  -m diamond \
    --data_dir /lustre/home/acct-ioozy/ioozy/software/eggnog-mapper/data --cpu 40
#TPM
sample_list=(Cadiz_Origin Cadiz_8MPa Cadiz_15MPa Cadiz_30MPa Cadiz_8MPaII Cadiz_15MPaII Cadiz_30MPaII)
for sample in ${sample_list[@]}
do
    gff='Cadiz_merge_03_gene_prodigal.gff'
    bam=${sample}_04_map_bbmap.bam
    featureCounts -a ${gff} -o 06_${sample}.count -t CDS -g ID -T 39 -p ${bam}
    04_tpm_featureCount.py 06_${sample}.count 07_${sample}.tpm
done
```

#### 5. MAG

```bash
ctg2mag=''
methods=''
sensitivity=(60 75 90)
specific=(60 75 90)
for p in ${sensitivity[@]}
do 
    for s in ${specific[@]}
    do
        mkdir metabat2_${p}_${s}
        metabat2 -i 03_spades.fa -a Cadiz_merge_05_depth_bbmap.metabat \
            -o metabat2_${p}_${s}/${p}_${s} --seed 599 \
            --maxP ${s} --minS ${p}
        21_MAGs_scaffold_list.py metabat2_${p}_${s} fa metabat2_${p}_${s}.txt
        ctg2mag="${ctg2mag},metabat2_${p}_${s}.txt"
        methods="${methods},${p}_${s}"
    done
done

cut -f 1,4 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_Origin.list
cut -f 1,6 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_8MPa.list
cut -f 1,8 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_15MPa.list
cut -f 1,10 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_30MPa.list
cut -f 1,12 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_8MPaII.list
cut -f 1,14 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_15MPaII.list
cut -f 1,16 Cadiz_merge_05_depth_bbmap.metabat > maxbin2_Cadiz_30MPaII.list
ls *list > maxbin2.abund

mkdir maxbin2_107
run_MaxBin.pl -contig 03_spades.fa \
    -markerset 107 -thread 40 \
    -out maxbin2_107/maxbin_107 \
    -abund_list maxbin2.abund 
21_MAGs_scaffold_list.py maxbin2_107 fasta maxbin2_107.txt
ctg2mag="${ctg2mag},maxbin2_107.txt"
methods="${methods},axbin2_107"

mkdir maxbin2_40 
run_MaxBin.pl -contig 03_spades.fa \
    -markerset 40 -thread 40 \
    -out maxbin2_40/maxbin_40 \
    -abund_list maxbin2.abund
21_MAGs_scaffold_list.py maxbin2_40 fasta maxbin2_40.txt
ctg2mag="${ctg2mag},maxbin2_40.txt"
methods="${methods},maxbin2_40"

cut_up_fasta.py 03_spades.fa -c 10000 -o 0 \
    --merge_last -b concoct_03_10K.bed > concoct_03_10K.fa
concoct_coverage_table.py concoct_03_10K.bed ../04_map_bbmap/*.bam > concoct_03_10K.tsv
concoct --coverage_file concoct_03_10K.tsv \
    --composition_file concoct_03_10K.fa \
    -t 40 -r 150 -b concoct -s 599 --no_original_data
merge_cutup_clustering.py concoct_clustering_gt1000.csv > concoct_clustering_gt1000_merge.csv
mkdir concoct_bins
extract_fasta_bins.py --output_path concoct_bins/ 03_spades.fa concoct_clustering_gt1000_merge.csv
21_MAGs_scaffold_list.py concoct_bins fa concoct.txt
ctg2mag="${ctg2mag},concoct.txt"
methods="${methods},concoct"

mkdir 01_DAS
DAS_Tool -i ${ctg2mag} -l ${methods} \
    -c 03_spades.fa -o 01_DAS/Cadiz  \
    --write_bins 1 -t 40 --search_engine diamond --score_threshold 0
    
checkm lineage_wf -x fa -t 40 ./01_DAS/Cadiz_DASTool_bins/ ./02_checkm 
##after correction 
gtdbtk classify_wf --cpus 39 --genome_dir ./00_MAGs -x fa --out_dir ./08_GTDBtk_r95
```

