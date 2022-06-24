# RNA-Seq

In this study, performed by Holzer et al. (2016), 
- (1) total RNA from a human HuH7 cell line and a fruit bat cell line (R06E-J; Rousettus aegyptiacs) infected with either the Ebola or Marburg virus (EBOV, MARV) was harvested 3, 7, and 23 h postinfection, depleted of ribosomal RNA and sequenced on an Illumina HiSeq2500. The bat RNA was further pooled and additionally sequenced on an Illumina MiSeq system. Initial quality control and trimming of the raw data were conducted with FastQC (Andrews, 2010) and PRINSEQ (Schmieder and Edwards, 2011). 
- (2) For bat RNA, a de novo transcriptome assembly was constructed by combining MiSeq and HiSeq data using Velvet/Oases (Schulz et al., 2012; Zerbino and Birney, 2008), ABySS/Trans-ABySS (Birol et al., 2009; Simpson et al., 2009), SOAPdenovo-Trans (Luo et al., 2012), Trinity (Grabherr et al., 2011), and Mira (Chevreux et al., 2004) with default parameters and multiple k-mer values, if possible. 
- (3) The mapping of the RNA-Seq short-reads was performed for Mock-, EBOV-, and MARV-treated cells onto human/bat genomes and the bat transcriptome with Segemehl (Hoffmann et al., 2014) and TopHat (Kim et al., 2013). 
- (4) A differential gene expression analysis was performed by counting uniquely mapped reads with HTSeq-count (Anders et al., 2015) and applying a DESeq (Love et al., 2014) analysis in R. The results were further used for clustering and scatter/group plot analyzes.
- (5) A homology search in bats was performed for all significantly differential expressed genes from (4) and for the genes assumed to be involved in the
response to infection based on an enriched pathway analysis and the literature. The Rousettus aegyptiacus genome and coding sequences from Pteropus vampyrus, a closely related bat species, were used to validate but also to detect homologous sequences in the bat transcriptome. Detected homologs were employed for the differential gene expression analysis. 
- (6) One huge advantage of this comprehensive study was the manual inspection of ~7.5 % of the human genes. Each candidate gene was manually investigated in the IGV (Thorvaldsdóttir et al., 2013) and UCSC (Dreszer et al., 2012) browsers for the human and bat samples from all time points. Single-nucleotide modifications (differential SNPs, posttranscriptional modifications), intronic transcripts and regulators, alternative splicing and isoforms, as well as upstream and downstream transcript characteristics were described.

## 0, create environment
conda env create -n rnaseq -c conda-forge -bioconda -defaults fastqc trim-galore star=2.6.1d

## 1(optional), mapping on the p602, p604, and p605 (MCPyV)
```sh
./damian_extended/main.rb --host human3 --type dna --ref_db MCPyV.fa --external_ref MCPyV.fa --external_groupref MCPyV.fa --min_contiglength 200  -u trimmed/V_8_4_2_p602_d8_DonorI.fastq.gz  --sample p602_on_MCPyV   --no_anno_groupref   --threads 15 --force

for rna_sample in p602_d8_DonorI p602_d8_DonorII p604_d8_DonorI p604_d8_DonorII p605_d8_DonorI p605_d8_DonorII; do
for rna_sample in p602_d8_DonorI p602_d8_DonorII; do
bwa mem -M -t 14 MCPyV.fa trimmed/${rna_sample}.fastq.gz > ${rna_sample}_on_MCPyV.sam
samtools view -Sb ${rna_sample}_on_MCPyV.sam > ${rna_sample}_on_MCPyV.bam
samtools flagstat ${rna_sample}_on_MCPyV.bam 
rm *.sam
done
```


## 2, raw data (from genbank or with own data)
```sh
# from genbank (see downloading_fastq_GEO.pdf)
cat SRR_Acc_List.txt | xargs -n 1 bash get_SRR_data.sh
gzip *.fastq

# with own data: RNASeq batch 1 --> good_quality
#RNASeq V.8.0
ln -s ../Raw_Data2/200420_NB501882_0195_AHJCNWBGXF/nf207/8_0_mock_NHDF_Donor_I_P9_S11_R1_001.fastq.gz  V_8_0_mock_DonorI.fastq.gz
ln -s ../Raw_Data2/200420_NB501882_0195_AHJCNWBGXF/nf208/8_0_mock_NHDF_Donor_II_P9_S12_R1_001.fastq.gz V_8_0_mock_DonorII.fastq.gz
#RNASeq V.8.1.5
ln -s ../Raw_Data2/190927_NB501882_0150_AHC2MHBGXC_single/nf135/8_1_5_p601_d3_S1_R1_001.fastq.gz V_8_1_5_p601_d3_DonorII.fastq.gz
ln -s ../Raw_Data2/190927_NB501882_0150_AHC2MHBGXC_single/nf136/8_1_5_p604_d3_S2_R1_001.fastq.gz V_8_1_5_p604_d3_DonorII.fastq.gz
ln -s ../Raw_Data2/190927_NB501882_0150_AHC2MHBGXC_single/nf137/8_1_5_p601_d8_S3_R1_001.fastq.gz V_8_1_5_p601_d8_DonorII.fastq.gz
ln -s ../Raw_Data2/190927_NB501882_0150_AHC2MHBGXC_single/nf138/8_1_5_p604_d8_S4_R1_001.fastq.gz V_8_1_5_p604_d8_DonorII.fastq.gz
#RNASeq V.8.1.3 (should be V.8.1.6)
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf29/8_1_3_p601_d3_S1_R1_001.fastq.gz   V_8_1_6_p601_d3_DonorI.fastq.gz
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf30/8_1_3_p601_d8_S2_R1_001.fastq.gz   V_8_1_6_p601_d8_DonorI.fastq.gz
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf31/_8_1_3_p604_d3_S3_R1_001.fastq.gz  V_8_1_6_p604_d3_DonorI.fastq.gz
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf32/8_1_3_p604_d8_S4_R1_001.fastq.gz   V_8_1_6_p604_d8_DonorI.fastq.gz
#RNASeq V.8.2.3
ln -s ../Raw_Data2/191008_NB501882_0151_AHC235BGXC/nf139/8_2_3_p600_d3_S20_R1_001.fastq.gz V_8_2_3_p600_d3_DonorII.fastq.gz
ln -s ../Raw_Data2/191008_NB501882_0151_AHC235BGXC/nf140/8_2_3_p605_d3_S21_R1_001.fastq.gz V_8_2_3_p605_d3_DonorII.fastq.gz
ln -s ../Raw_Data2/191008_NB501882_0151_AHC235BGXC/nf141/8_2_3_p600_d8_S22_R1_001.fastq.gz V_8_2_3_p600_d8_DonorII.fastq.gz
ln -s ../Raw_Data2/191008_NB501882_0151_AHC235BGXC/nf142/8_2_3_p605_d8_S23_R1_001.fastq.gz V_8_2_3_p605_d8_DonorII.fastq.gz
#RNASeq V.8.2.1 (should be V.8.2.4)
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf33/8_2_1_p600_d3_S5_R1_001.fastq.gz   V_8_2_4_p600_d3_DonorI.fastq.gz
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf34/8_2_1_p600_d8_S6_R1_001.fastq.gz   V_8_2_4_p600_d8_DonorI.fastq.gz
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf35/8_2_1_p605_d3_S7_R1_001.fastq.gz   V_8_2_4_p605_d3_DonorI.fastq.gz
ln -s ../Raw_Data2/190228_NB501882_0102_AHYJJ2BGX9/nf36/8_2_1_p605_d8_S8_R1_001.fastq.gz   V_8_2_4_p605_d8_DonorI.fastq.gz
#RNASeq V.8.4.1 (Note V_8_4_1 should be V_8_4_2)
ln -s ../Raw_Data2/200817_NB501882_0221_AHFMVLBGXG/nf298/V_8_4_1_p602_Donor_II_d8_S12_R1_001.fastq.gz V_8_4_2_p602_d8_DonorII.fastq.gz
ln -s ../Raw_Data2/200817_NB501882_0221_AHFMVLBGXG/nf297/V_8_4_2_p602_Donor_I_d8_S11_R1_001.fastq.gz  V_8_4_2_p602_d8_DonorI.fastq.gz
#RNASeq V.8.4.2 (Note V_8_4_2 should be V_8_4_1)
ln -s ../Raw_Data2/210302_NB501882_0250_AHLM2YBGXH/nf613/V842_NHDF_Donor_1_LT_d3_S19_R1_001.fastq.gz  V_8_4_1_p602_d3_DonorI.fastq.gz
ln -s ../Raw_Data2/210302_NB501882_0250_AHLM2YBGXH/nf614/V842_NHDF_Donor_2_LT_d3_S20_R1_001.fastq.gz  V_8_4_1_p602_d3_DonorII.fastq.gz
#RNASeq V.8.3.1 (Note that 8_2_3 should be V_8_3_1)
ln -s ../Raw_Data2/200420_NB501882_0195_AHJCNWBGXF/nf209/8_3_1_p600and601_d12_S13_R1_001.fastq.gz    V_8_3_1_p600and601_d12_DonorI.fastq.gz
ln -s ../Raw_Data2/200420_NB501882_0195_AHJCNWBGXF/nf210/8_2_3_p604and605_d12_S14_R1_001.fastq.gz    V_8_3_1_p604and605_d12_DonorI.fastq.gz
#RNASeq V.8.3.2
ln -s ../Raw_Data2/200817_NB501882_0221_AHFMVLBGXG/nf299/V_8_3_2_p600_601_Donor_II_d9_S13_R1_001.fastq.gz  V_8_3_2_p600and601_d9_DonorII.fastq.gz
ln -s ../Raw_Data2/200817_NB501882_0221_AHFMVLBGXG/nf300/V_8_3_2_p604_605_Donor_II_d9_S14_R1_001.fastq.gz  V_8_3_2_p604and605_d9_DonorII.fastq.gz

#---- to compare ----
P600 vs untreated / mock (donor 1+2)
P601 vs untreated / mock (donor 1+2)
P604 d3 vs p601 d3 (Donor 1+2)
P604 d8 vs p601 d8 (donor 1+2)
P605 d3 vs p600 d3 (donor 1+2)
P605 d8 vs p600 d8 (donor 1+2)
(These were already analysed)

What is new:
P602 d8 vs p600 d8 (donor 1+2)
P604+605 d9/12 vs p600+601 d9/12 (donor 1+2)
P602 d3 vs p600 d3 (donor 1+2)

./190228_NB501882_0102_AHYJJ2BGX9/nf32/8_1_3_p604_d8_S4_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf29/8_1_3_p601_d3_S1_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf31/_8_1_3_p604_d3_S3_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf36/8_2_1_p605_d8_S8_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf30/8_1_3_p601_d8_S2_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf34/8_2_1_p600_d8_S6_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf35/8_2_1_p605_d3_S7_R1_001.fastq.gz
./190228_NB501882_0102_AHYJJ2BGX9/nf33/8_2_1_p600_d3_S5_R1_001.fastq.gz

P601: mock (sT)
P604: sT
P600: mock (LT trunc.)
P605: LT trunc.

#--batch1 (190228_NB501882_0102_AHYJJ2BGX9)--
mv 8_1_3_p601_d3_S1_R1_001.fastq.gz mock_sT_d3.fastq.gz
mv 8_1_3_p601_d8_S2_R1_001.fastq.gz mock_sT_d8.fastq.gz
mv 8_1_3_p604_d3_S3_R1_001.fastq.gz sT_d3.fastq.gz
mv 8_1_3_p604_d8_S4_R1_001.fastq.gz sT_d8.fastq.gz
mv 8_2_1_p600_d3_S5_R1_001.fastq.gz mock_truncLT_d3.fastq.gz
mv 8_2_1_p600_d8_S6_R1_001.fastq.gz mock_truncLT_d8.fastq.gz
mv 8_2_1_p605_d3_S7_R1_001.fastq.gz truncLT_d3.fastq.gz
mv 8_2_1_p605_d8_S8_R1_001.fastq.gz truncLT_d8.fastq.gz
#--batch2 (190927_NB501882_0150_AHC2MHBGXC_single) --
mv 8_1_5_p601_d3_S1_R1_001.fastq.gz mock_sT_d3_r2.fastq.gz
mv 8_1_5_p601_d8_S3_R1_001.fastq.gz mock_sT_d8_r2.fastq.gz
mv 8_1_5_p604_d3_S2_R1_001.fastq.gz sT_d3_r2.fastq.gz
mv 8_1_5_p604_d8_S4_R1_001.fastq.gz sT_d8_r2.fastq.gz
#--batch3 (191008_NB501882_0151_AHC235BGXC) --
mv 8_2_3_p600_d3_S20_R1_001.fastq.gz mock_truncLT_d3_r2.fastq.gz
mv 8_2_3_p600_d8_S22_R1_001.fastq.gz mock_truncLT_d8_r2.fastq.gz
mv 8_2_3_p605_d3_S21_R1_001.fastq.gz truncLT_d3_r2.fastq.gz
mv 8_2_3_p605_d8_S23_R1_001.fastq.gz truncLT_d8_r2.fastq.gz
```


## 3, trimming (included in nextflow, not necessary) + nextflow
```sh
#http://genome.ucsc.edu/FAQ/FAQformat.html#format4
--0--
#mkdir trimmed; 
#cd trimmed;
#for sample_id in mock_sT_d3 mock_sT_d8 sT_d3 sT_d8 mock_truncLT_d3 mock_truncLT_d8 truncLT_d3 truncLT_d8 truncLT_d8_r2 mock_truncLT_d3_r2 truncLT_d3_r2 mock_truncLT_d8_r2 sT_d8_r2 mock_sT_d3_r2 sT_d3_r2 mock_sT_d8_r2    IMR90_GFP_0_1 IMR90_GFP_0_2 IMR90_GFP_0_3 IMR90_GFP_8_1 IMR90_GFP_8_2 IMR90_GFP_8_3 IMR90_GFP_16_1 IMR90_GFP_16_2 IMR90_GFP_16_3 IMR90_GFP_24_1 IMR90_GFP_24_2 IMR90_GFP_24_3 IMR90_GFP_32_1 IMR90_GFP_32_2 IMR90_GFP_32_3 IMR90_GFP_40_1 IMR90_GFP_40_2 IMR90_GFP_40_3 IMR90_GFP_48_1 IMR90_GFP_48_2 IMR90_GFP_48_3 IMR90_GFP_56_1 IMR90_GFP_56_2 IMR90_GFP_56_3 IMR90_GFP_64_1 IMR90_GFP_64_2 IMR90_GFP_64_3 IMR90_GFP_72_1 IMR90_GFP_72_2 IMR90_GFP_72_3 IMR90_GFP_80_1 IMR90_GFP_80_2 IMR90_GFP_80_3 IMR90_GFP_88_1 IMR90_GFP_88_2 IMR90_GFP_88_3 IMR90_GFP_96_1 IMR90_GFP_96_2 IMR90_GFP_96_3 IMR90_MCPyV_ST_0_1 IMR90_MCPyV_ST_0_2 IMR90_MCPyV_ST_0_3 IMR90_MCPyV_ST_8_1 IMR90_MCPyV_ST_8_2 IMR90_MCPyV_ST_8_3 IMR90_MCPyV_ST_16_1 IMR90_MCPyV_ST_16_2 IMR90_MCPyV_ST_16_3 IMR90_MCPyV_ST_24_1 IMR90_MCPyV_ST_24_2 IMR90_MCPyV_ST_24_3 IMR90_MCPyV_ST_32_1 IMR90_MCPyV_ST_32_2 IMR90_MCPyV_ST_32_3 IMR90_MCPyV_ST_40_1 IMR90_MCPyV_ST_40_2 IMR90_MCPyV_ST_40_3 IMR90_MCPyV_ST_48_1 IMR90_MCPyV_ST_48_2 IMR90_MCPyV_ST_48_3 IMR90_MCPyV_ST_56_1 IMR90_MCPyV_ST_56_2 IMR90_MCPyV_ST_56_3 IMR90_MCPyV_ST_64_1 IMR90_MCPyV_ST_64_2 IMR90_MCPyV_ST_64_3 IMR90_MCPyV_ST_72_1 IMR90_MCPyV_ST_72_2 IMR90_MCPyV_ST_72_3 IMR90_MCPyV_ST_80_1 IMR90_MCPyV_ST_80_2 IMR90_MCPyV_ST_80_3 IMR90_MCPyV_ST_88_1 IMR90_MCPyV_ST_88_2 IMR90_MCPyV_ST_88_3 IMR90_MCPyV_ST_96_1 IMR90_MCPyV_ST_96_2 IMR90_MCPyV_ST_96_3; do \
#        java -jar /home/jhuang/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 16 ../Raw_Data/${sample_id}.fastq.gz ${sample_id}.fastq.gz ILLUMINACLIP:/home/jhuang/Tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 AVGQUAL:20; \
#done 2>trimmomatic.log

# see read_biotype_assignment: https://github.com/ewels/ngi_visualizations/tree/master/ngi_visualizations/biotypes
#Running twice, once without biotype using UCSC  --> get the results which comparable to ChIPSeq-data, since they both using UCSC-gtf -->1
               once with biotype using ENSEMBLE reference --> get images of biotypes and as controls -->2

#--1--
## since in UCSC-gtf, there is no gene_biotype, we use gene_id in --fcGroupFeaturesType; NOTE that bed-file and HISAT-index are generated on site!
#nextflow run rnaseq --reads '/home/jhuang/DATA/Data_Denise_RNASeq/Raw_Data/*.fastq.gz' --fasta "/ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"  #--singleEnd -profile standard --aligner hisat2 --fcGroupFeaturesType gene_id --outdir results_nfcore -w work_nfcore -resume
## when using already generated genes.bed and HISAT2Index
#--hisat2_index "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/HISAT2Index/" --bed12 "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"
#
#--2--
## default type of --fcGroupFeaturesType is 'gene_biotype' if using Ensembl-gtf
#nextflow run rnaseq --reads '/home/jhuang/DATA/Data_Denise_RNASeq/Raw_Data/*.fastq.gz' --fasta "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa" --gtf "/media/jhuang/Titisee/ref/#Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf" --singleEnd -profile standard --aligner hisat2      -resume



#### prefer to use STARIndex, since HISAT2 needs memory of 120GB, otherwise it does not work for splicing mapping.
--1--
# since in UCSC-gtf, there is no gene_biotype, we use gene_id in --fcGroupFeaturesType --> NOTE that bed-file and HISAT-index are generated on site!
#DEBUG: An example of attributes included in your GTF annotation is 'gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16932";'
#---- SUCCESSFUL using STAR ----
ln -s /home/jhuang/Tools/rnaseq rnaseq
nextflow run rnaseq --reads "/home/jhuang/DATA/Data_Anastasia_RNASeq/Raw_Data/*.fastq.gz" --fasta "/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf" --star_index "/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/" --bed12 "/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"  --singleEnd -profile standard --aligner star --fcGroupFeaturesType gene_biotype --skip_genebody_coverage -resume

nextflow run rnaseq --reads "/home/jhuang/DATA/Data_Ute_RNASeq/Raw_Data/*.fastq.gz" --fasta "/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf" --star_index "/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/" --bed12 "/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"  --singleEnd -profile standard --aligner star --fcGroupFeaturesType gene_biotype --skip_genebody_coverage -resume

nextflow run rnaseq --reads "/home/jhuang/DATA/Data_Denise_RNASeq_plus_GSE79958/Raw_Data/*.fastq.gz" --fasta "/ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf" --star_index "/ref/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/" --bed12 "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed"  --singleEnd -profile standard --aligner star --fcGroupFeaturesType gene_id  -resume

#set pbs as the process executors.
nextflow run rnaseq --reads "/home/jhuang/DATA/Data_Denise_RNASeq/Raw_Data/*.fastq.gz" --fasta "/ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf" --star_index "/ref/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/" --bed12 "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed"  --singleEnd -profile standard --aligner star --fcGroupFeaturesType gene_id  -resume

#when using already generated genes.bed and *Index
--hisat2_index "/ref/Homo_sapiens/UCSC/hg19/Sequence/HISAT2Index/" (todo)   --bed12 "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed"
--star_index "/ref/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/" --bed12 "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed"

--2--
#default type of --fcGroupFeaturesType is 'gene_biotype' if using Ensembl-gtf
nextflow run rnaseq --reads '/home/jhuang/DATA/Data_Denise_RNASeq/Raw_Data/*.fastq.gz' --fasta "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa" --gtf "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf" --singleEnd -profile standard --aligner star   --outdir results_gr_gene_biotype -w work_gr_gene_biotype   -resume
#when using already generated genes.bed and *Index
--hisat2_index "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/HISAT2Index/" --bed12 "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"
--star_index "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/" (todo)    --bed12 "/media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"



#-- BUGs while preparing gtf and bed --
#under ngi_chipseq_ac2 using HISAT2
#--star_index "/ref/Homo_sapiens/UCSC/hg19/Sequence/HISAT2Index/" --> DEBUG3
nextflow run rnaseq --reads '/home/jhuang/DATA/Data_Denise_RNASeq/Raw_Data/*.fastq.gz'  --fasta "/ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf" --bed12 "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed" --singleEnd -profile standard --aligner hisat2  -resume
#under ngi_chipseq_ac2 using STAR
#nextflow run rnaseq --reads '/home/jhuang/DATA/Data_Denise_RNASeq/Raw_Data/*.fastq.gz'  --star_index "/ref/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/" --fasta "/ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" --gtf "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"  --bed12 "/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed" --singleEnd -profile standard     -resume
--downloadGTF 
https://www.ensembl.org/info/data/ftp/index.html
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
http://genome.ucsc.edu/cgi-bin/hgTables
#DEBUG1
conda install -c bioconda bioconductor-rsubread
conda install -c bioconda subread
conda install -c bioconda stringtie
conda install -c bioconda hisat2
conda update hisat2
conda install -c bioconda preseq
# installed Rsubread under ~/anaconda3/envs/ngi_chipseq_ac2/lib/R/library
BiocManager::install(c("Rsubread", "dupRadar", "limma", "lattice", "locfit", "edgeR", "chron", "data.table", "gtools", "gdata", "bitops", "caTools", "gplots", "markdown"))
#DEBUG2
 ERROR: failed to find the gene identifier attribute in the 9th column of the provided GTF file.
  The specified gene identifier attribute is 'gene_biotype'
  An example of attributes included in your GTF annotation is 'gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16932";'
  The program has to terminate.
featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
featureCounts -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n 7 > tmp_file
cat $biotypes_header tmp_file >> ${bam_featurecounts.baseName}_biotype_counts_mqc.txt
#needs the gtf of ENSEMBLE, since it includes "gene_biotype" and "gene_id"
# ENSEMBLE (874M) 9th colum
1 processed_transcript                transcript  11869 14409 . + . gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_sourc e "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
# UCSC (136M) 9th column
gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS8568";
diff /ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
cut -f1-2 /ref/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai | sort > hg19_UCSC.fa_length
cut -f1-2 /media/jhuang/Titisee/ref/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa.fai | sort > hg19_Ensemble.fa_length
cut -f2 hg19_UCSC.fa_length > hg19_UCSC.fa_length_
cut -f2 hg19_Ensemble.fa_length > hg19_Ensemble.fa_length_



#CONSOLE
cut -f2- merged_gene_counts.txt > merged_gene_counts_2.txt
```


## 3.5, rerun multiqc
```sh
mv MultiQC .. 
mv pipeline_info .. 
mv STAR ..
mv rseqc/read_distribution ..
mv markDuplicates/metrics ..
#Note that sequence counts occurs twice in the results, featureCounts (16,9m/68.0%=24,85m), Biotype_Counts(17 m) and FastQC/Sequence_Counts (16,03m/74,8%=21,4m after quality control and filtering).
#Note that Dups only uses the FastQC-results
cp ~/DATA/results_description.html ./

ln -s ~/Tools/rnaseq/assets/multiqc_config.yaml multiqc_config.yaml
multiqc -f --config multiqc_config.yaml . 2>&1
rm multiqc_config.yaml
```


## 4, load and clean data, and construct DESeqDataSet
```sh
library("AnnotationDbi")
library("clusterProfiler")
library("ReactomePA")
library("org.Hs.eg.db")
library(DESeq2)
library(gplots)
setwd("/mnt/Seagate_Corona/Data_Denise_RNASeq/results/featureCounts")

####-- dataset 16 isolates: CORRECTION!! --

# --- IMPORTANT: to reminder that sequencing in a batch if possible! ---
d.raw<- read.delim2("merged_gene_counts_2.txt",sep="\t", header=TRUE, row.names=1)
 [1] V_8_0_mock_DonorIAligned.sortedByCoord.out.bam            
 [2] V_8_1_6_p601_d8_DonorIAligned.sortedByCoord.out.bam       
 [3] V_8_1_5_p604_d3_DonorIIAligned.sortedByCoord.out.bam      
 [4] V_8_1_6_p604_d3_DonorIAligned.sortedByCoord.out.bam       
 [5] V_8_2_3_p605_d8_DonorIIAligned.sortedByCoord.out.bam      
 [6] V_8_0_mock_DonorIIAligned.sortedByCoord.out.bam           
 [7] V_8_1_5_p601_d8_DonorIIAligned.sortedByCoord.out.bam      
 [8] V_8_2_3_p605_d3_DonorIIAligned.sortedByCoord.out.bam      
 [9] V_8_1_6_p604_d8_DonorIAligned.sortedByCoord.out.bam       
[10] V_8_2_3_p600_d8_DonorIIAligned.sortedByCoord.out.bam      
[11] V_8_2_4_p600_d8_DonorIAligned.sortedByCoord.out.bam       
[12] V_8_2_4_p600_d3_DonorIAligned.sortedByCoord.out.bam       
[13] V_8_1_5_p604_d8_DonorIIAligned.sortedByCoord.out.bam      
[14] V_8_1_5_p601_d3_DonorIIAligned.sortedByCoord.out.bam      
[15] V_8_1_6_p601_d3_DonorIAligned.sortedByCoord.out.bam       
[16] V_8_2_3_p600_d3_DonorIIAligned.sortedByCoord.out.bam      
[17] V_8_2_4_p605_d3_DonorIAligned.sortedByCoord.out.bam       
[18] V_8_4_1_p602_d8_DonorIAligned.sortedByCoord.out.bam       
[19] V_8_3_2_p604and605_d9_DonorIIAligned.sortedByCoord.out.bam
[20] V_8_4_2_p602_d3_DonorIAligned.sortedByCoord.out.bam       
[21] V_8_4_2_p602_d3_DonorIIAligned.sortedByCoord.out.bam      
[22] V_8_2_4_p605_d8_DonorIAligned.sortedByCoord.out.bam       
[23] V_8_3_1_p600and601_d12_DonorIAligned.sortedByCoord.out.bam
[24] V_8_4_1_p602_d8_DonorIIAligned.sortedByCoord.out.bam      
[25] V_8_3_1_p604and605_d12_DonorIAligned.sortedByCoord.out.bam
[26] V_8_3_2_p600and601_d9_DonorIIAligned.sortedByCoord.out.bam

colnames(d.raw)<- c("V_8_0_mock_DonorI","V_8_1_6_p601_d8_DonorI","V_8_1_5_p604_d3_DonorII","V_8_1_6_p604_d3_DonorI","V_8_2_3_p605_d8_DonorII","V_8_0_mock_DonorII","V_8_1_5_p601_d8_DonorII","V_8_2_3_p605_d3_DonorII","V_8_1_6_p604_d8_DonorI","V_8_2_3_p600_d8_DonorII","V_8_2_4_p600_d8_DonorI","V_8_2_4_p600_d3_DonorI","V_8_1_5_p604_d8_DonorII","V_8_1_5_p601_d3_DonorII","V_8_1_6_p601_d3_DonorI","V_8_2_3_p600_d3_DonorII","V_8_2_4_p605_d3_DonorI","V_8_4_1_p602_d8_DonorI","V_8_3_2_p604and605_d9_DonorII","V_8_4_2_p602_d3_DonorI","V_8_4_2_p602_d3_DonorII","V_8_2_4_p605_d8_DonorI","V_8_3_1_p600and601_d12_DonorI","V_8_4_1_p602_d8_DonorII","V_8_3_1_p604and605_d12_DonorI","V_8_3_2_p600and601_d9_DonorII")  #26364

col_order <- c("V_8_0_mock_DonorI","V_8_0_mock_DonorII","V_8_1_5_p601_d3_DonorII", "V_8_1_5_p604_d3_DonorII", "V_8_1_5_p601_d8_DonorII","V_8_1_5_p604_d8_DonorII",   "V_8_1_6_p601_d3_DonorI","V_8_1_6_p604_d3_DonorI","V_8_1_6_p601_d8_DonorI","V_8_1_6_p604_d8_DonorI",  "V_8_2_3_p600_d3_DonorII","V_8_2_3_p605_d3_DonorII","V_8_2_3_p600_d8_DonorII", "V_8_2_3_p605_d8_DonorII",  "V_8_2_4_p600_d3_DonorI","V_8_2_4_p605_d3_DonorI","V_8_2_4_p600_d8_DonorI","V_8_2_4_p605_d8_DonorI",  "V_8_4_1_p602_d8_DonorII","V_8_4_1_p602_d8_DonorI",  "V_8_3_1_p600and601_d12_DonorI", "V_8_3_1_p604and605_d12_DonorI","V_8_3_2_p600and601_d9_DonorII","V_8_3_2_p604and605_d9_DonorII",    "V_8_4_2_p602_d3_DonorI","V_8_4_2_p602_d3_DonorII")
reordered.raw <- d.raw[,col_order]
#reordered.raw <- subset(d.raw, select=col_order)

## OPTION1: sT
#d.raw_sT <- subset(d.raw, select=c("mock_sT_d3_r2","mock_sT_d8","sT_d3","sT_d8_r2","mock_sT_d3","sT_d3_r2","mock_sT_d8_r2","sT_d8"))
#d <- d.raw_sT[rowSums(d.raw_sT>3)>2,]
#replicates = as.factor(c("mock_sT_d3","mock_sT_d8","sT_d3","sT_d8","mock_sT_d3","sT_d3","mock_sT_d8","sT_d8"))
#batch = as.factor(c("r2","r1","r1","r2","r1","r2","r2","r1"))
#ids = as.factor(c("mock_sT_d3_r2","mock_sT_d8","sT_d3","sT_d8_r2","mock_sT_d3","sT_d3_r2","mock_sT_d8_r2","sT_d8"))
#cData = data.frame(row.names=colnames(d), replicates=replicates, batch=batch, ids=ids)
#dds<-DESeqDataSetFromMatrix(countData=d, colData=cData, design=~batch+replicates)

## OPTION2: truncLT
#d.raw_truncLT <- subset(d.raw, #select=c("truncLT_d8_r2","truncLT_d3","mock_truncLT_d8_r2","mock_truncLT_d3_r2","mock_truncLT_d8","mock_truncLT_d3","truncLT_d3_r2","truncLT_d8"))
#d <- d.raw_truncLT[rowSums(d.raw_truncLT>3)>2,]
#replicates = as.factor(c("truncLT_d8","truncLT_d3","mock_truncLT_d8","mock_truncLT_d3","mock_truncLT_d8","mock_truncLT_d3","truncLT_d3","truncLT_d8"))
#batch = as.factor(c("r2","r1","r2","r2","r1","r1","r2","r1"))
#ids = as.factor(c("truncLT_d8_r2","truncLT_d3","mock_truncLT_d8_r2","mock_truncLT_d3_r2","mock_truncLT_d8","mock_truncLT_d3","truncLT_d3_r2","truncLT_d8"))
#cData = data.frame(row.names=colnames(d), replicates=replicates, batch=batch, ids=ids)
#dds<-DESeqDataSetFromMatrix(countData=d, colData=cData, design=~batch+replicates)

## OPTION3: IPA sT
#d.raw_sT <- subset(d.raw, select=c("mock_sT_d3_r2","mock_sT_d8","sT_d3","sT_d8_r2","mock_sT_d3","sT_d3_r2","mock_sT_d8_r2","sT_d8"))
#replicates = as.factor(c("mock_sT_d3","mock_sT_d8","sT_d3","sT_d8","mock_sT_d3","sT_d3","mock_sT_d8","sT_d8"))
#batch = as.factor(c("r2","r1","r1","r2","r1","r2","r2","r1"))
#ids = as.factor(c("mock_sT_d3_r2","mock_sT_d8","sT_d3","sT_d8_r2","mock_sT_d3","sT_d3_r2","mock_sT_d8_r2","sT_d8"))
#cData = data.frame(row.names=colnames(d.raw_sT), replicates=replicates, batch=batch, ids=ids)
#dds<-DESeqDataSetFromMatrix(countData=d.raw_sT, colData=cData, design=~batch+replicates)

## OPTION4: IPA truncLT
#d.raw_truncLT <- subset(d.raw, #select=c("truncLT_d8_r2","truncLT_d3","mock_truncLT_d8_r2","mock_truncLT_d3_r2","mock_truncLT_d8","mock_truncLT_d3","truncLT_d3_r2","truncLT_d8"))
#replicates = as.factor(c("truncLT_d8","truncLT_d3","mock_truncLT_d8","mock_truncLT_d3","mock_truncLT_d8","mock_truncLT_d3","truncLT_d3","truncLT_d8"))
#batch = as.factor(c("r2","r1","r2","r2","r1","r1","r2","r1"))
#ids = as.factor(c("truncLT_d8_r2","truncLT_d3","mock_truncLT_d8_r2","mock_truncLT_d3_r2","mock_truncLT_d8","mock_truncLT_d3","truncLT_d3_r2","truncLT_d8"))
#cData = data.frame(row.names=colnames(d.raw_truncLT), replicates=replicates, batch=batch, ids=ids)
#dds<-DESeqDataSetFromMatrix(countData=d.raw_truncLT, colData=cData, design=~batch+replicates)

## OPTION5: sT vs untreated or LTtr vs untreated
#d.raw_mock_untreated <- subset(d.raw, select=c("mock_sT_d3_r2","mock_sT_d8_r1","mock_sT_d8_r2","mock_truncLT_d8_r2","mock_truncLT_d3_r1","mock_truncLT_d8_r1","mock_truncLT_d3_r2","mock_sT_d3_r1","untreated_r1","untreated_r2"))
#replicates = as.factor(c("sT_mock_d3","sT_mock_d8","sT_mock_d8","LTtr_mock_d8","LTtr_mock_d3","LTtr_mock_d8","LTtr_mock_d3","sT_mock_d3","untreated","untreated"))
#batch = as.factor(c("r2","r1","r2","r2","r1","r1","r2","r1","r3","r3"))
#ids = as.factor(c("sT_mock_d3_r2","sT_mock_d8_r1","sT_mock_d8_r2","LTtr_mock_d8_r2","LTtr_mock_d3_r1","LTtr_mock_d8_r1","LTtr_mock_d3_r2","sT_mock_d3_r1","untreated_r1","untreated_r2"))
#cData = data.frame(row.names=colnames(d.raw_mock_untreated), replicates=replicates, batch=batch, ids=ids)
#dds<-DESeqDataSetFromMatrix(countData=d.raw_mock_untreated, colData=cData, design=~batch+replicates)

# OPTION6: sT vs untreated or LTtr vs untreated
#P602 d8 vs p600 d8 (donor 1+2)
#P604+605 d9/12 vs p600+601 d9/12 (donor 1+2)
d.raw_p602_p600_mixed <- subset(reordered.raw, select=c("V_8_4_1_p602_d8_DonorII","V_8_4_2_p602_d8_DonorI","V_8_2_3_p600_d8_DonorII", "V_8_2_1_p600_d8_DonorI",   "V_8_3_2_p604and605_d9_DonorII","V_8_3_1_p604and605_d12_DonorI",  "V_8_3_2_p600and601_d9_DonorII", "V_8_3_1_p600and601_d12_DonorI"))
replicates = as.factor(c("p602_d8","p602_d8","p600_d8","p600_d8", "p604and605_d9ord12","p604and605_d9ord12","p600and601_d9ord12","p600and601_d9ord12"))
batch = as.factor(c("200817","200817","191008","190228","200817","200420","100817","200420"))
ids = as.factor(c("p602_d8_DII","p602_d8_DI","p600_d8_DII","p600_d8_DI", "p604and605_d9ord12","p604and605_d9ord12", "p600and601_d9ord12","p600and601_d9ord12"))
cData = data.frame(row.names=colnames(d.raw_p602_p600_mixed), replicates=replicates, batch=batch, ids=ids)
dds<-DESeqDataSetFromMatrix(countData=d.raw_p602_p600_mixed, colData=cData, design=~batch+replicates)

# OPTION7
#replicates = as.factor(c("untreated","untreated","p601_d3","p604_d3",     "p601_d8","p604_d8","p601_d3","p604_d3",     "p601_d8","p604_d8","p600_d3","p605_d3",   "p600_d8", "p605_d8","p600_d3","p605_d3",       "p600_d8","p605_d8","p602_d8","p602_d8",      "p600and601_d12", "p604and605_d12","p600and601_d9","p604and605_d9"))

replicates = as.factor(c("untreated","untreated","p601_d3","p604_d3",     "p601_d8","p604_d8","p601_d3","p604_d3",     "p601_d8","p604_d8","p600_d3","p605_d3",   "p600_d8", "p605_d8","p600_d3","p605_d3",       "p600_d8","p605_d8","p602_d8","p602_d8",      "p600and601_d9or12", "p604and605_d9or12","p600and601_d9or12","p604and605_d9or12"))
batch = as.factor(c("200420", "200420", "190927", "190927",    "190927", "190927", "190228", "190228",    "190228", "190228", "191008", "191008",    "191008", "191008", "190228", "190228",     "190228", "190228", "200817", "200817",       "200420", "200420", "200817", "200817"))
ids = as.factor(c("untreated_DonorI","untreated_DonorII", "p601_d3_DonorII","p604_d3_DonorII", "p601_d8_DonorII","p604_d8_DonorII",   "p601_d3_DonorI","p604_d3_DonorI","p601_d8_DonorI","p604_d8_DonorI",  "p600_d3_DonorII","p605_d3_DonorII","p600_d8_DonorII", "p605_d8_DonorII",  "p600_d3_DonorI","p605_d3_DonorI","p600_d8_DonorI","p605_d8_DonorI",  "p602_d8_DonorII","p602_d8_DonorI",  "p600and601_d12_DonorI", "p604and605_d12_DonorI","p600and601_d9_DonorII","p604and605_d9_DonorII"))
donor = as.factor(c("DonorI","DonorII", "DonorII","DonorII", "DonorII","DonorII",   "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorII","DonorII", "DonorII",  "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorI",  "DonorI", "DonorI","DonorII","DonorII"))
#Note that we need reordered.raw
cData = data.frame(row.names=colnames(reordered.raw), replicates=replicates, donor=donor, batch=batch, ids=ids)
dds<-DESeqDataSetFromMatrix(countData=reordered.raw, colData=cData, design=~donor+replicates)  #batch+
```

## 5, fake or non-fake replicates
```sh
#"p602_d3_DonorI","p602_d3_DonorII"
#, "210302", "210302"
# ---- ACTIVATED for fake replicates, y are the fake-replicate of x ----
# fake replicates
fake_reordered.raw <- reordered.raw+1
total <- merge(reordered.raw,fake_reordered.raw, by=0)
total2 <- total[,-1]
rownames(total2) <- total[,1]

replicates = as.factor(c("untreated_DonorI","untreated_DonorII", "p601_d3_DonorII","p604_d3_DonorII", "p601_d8_DonorII","p604_d8_DonorII",   "p601_d3_DonorI","p604_d3_DonorI","p601_d8_DonorI","p604_d8_DonorI",  "p600_d3_DonorII","p605_d3_DonorII","p600_d8_DonorII", "p605_d8_DonorII",  "p600_d3_DonorI","p605_d3_DonorI","p600_d8_DonorI","p605_d8_DonorI",  "p602_d8_DonorII","p602_d8_DonorI",  "p600and601_d12_DonorI", "p604and605_d12_DonorI","p600and601_d9_DonorII","p604and605_d9_DonorII",        "untreated_DonorI","untreated_DonorII", "p601_d3_DonorII","p604_d3_DonorII", "p601_d8_DonorII","p604_d8_DonorII",   "p601_d3_DonorI","p604_d3_DonorI","p601_d8_DonorI","p604_d8_DonorI",  "p600_d3_DonorII","p605_d3_DonorII","p600_d8_DonorII", "p605_d8_DonorII",  "p600_d3_DonorI","p605_d3_DonorI","p600_d8_DonorI","p605_d8_DonorI",  "p602_d8_DonorII","p602_d8_DonorI",  "p600and601_d12_DonorI", "p604and605_d12_DonorI","p600and601_d9_DonorII","p604and605_d9_DonorII"))

batch = as.factor(c("200420", "200420", "190927", "190927",    "190927", "190927", "190228", "190228",    "190228", "190228", "191008", "191008",    "191008", "191008", "190228", "190228",     "190228", "190228", "200817", "200817",       "200420", "200420", "200817", "200817",       "200420", "200420", "190927", "190927",    "190927", "190927", "190228", "190228",    "190228", "190228", "191008", "191008",    "191008", "191008", "190228", "190228",     "190228", "190228", "200817", "200817",       "200420", "200420", "200817", "200817"))

ids = as.factor(c("untreated_DonorI.x","untreated_DonorII.x", "p601_d3_DonorII.x","p604_d3_DonorII.x", "p601_d8_DonorII.x","p604_d8_DonorII.x",   "p601_d3_DonorI.x","p604_d3_DonorI.x","p601_d8_DonorI.x","p604_d8_DonorI.x",  "p600_d3_DonorII.x","p605_d3_DonorII.x","p600_d8_DonorII.x", "p605_d8_DonorII.x",  "p600_d3_DonorI.x","p605_d3_DonorI.x","p600_d8_DonorI.x","p605_d8_DonorI.x",  "p602_d8_DonorII.x","p602_d8_DonorI.x",  "p600and601_d12_DonorI.x", "p604and605_d12_DonorI.x","p600and601_d9_DonorII.x","p604and605_d9_DonorII.x",        "untreated_DonorI.y","untreated_DonorII.y", "p601_d3_DonorII.y","p604_d3_DonorII.y", "p601_d8_DonorII.y","p604_d8_DonorII.y",   "p601_d3_DonorI.y","p604_d3_DonorI.y","p601_d8_DonorI.y","p604_d8_DonorI.y",  "p600_d3_DonorII.y","p605_d3_DonorII.y","p600_d8_DonorII.y", "p605_d8_DonorII.y",  "p600_d3_DonorI.y","p605_d3_DonorI.y","p600_d8_DonorI.y","p605_d8_DonorI.y",  "p602_d8_DonorII.y","p602_d8_DonorI.y",  "p600and601_d12_DonorI.y", "p604and605_d12_DonorI.y","p600and601_d9_DonorII.y","p604and605_d9_DonorII.y", ))

donor = as.factor(c("DonorI","DonorII", "DonorII","DonorII", "DonorII","DonorII",   "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorII","DonorII", "DonorII",  "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorI",  "DonorI", "DonorI","DonorII","DonorII",        "DonorI","DonorII", "DonorII","DonorII", "DonorII","DonorII",   "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorII","DonorII", "DonorII",  "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorI",  "DonorI", "DonorI","DonorII","DonorII"))

cData = data.frame(row.names=colnames(total2), replicates=replicates, donor=donor, batch=batch, ids=ids)
dds<-DESeqDataSetFromMatrix(countData=total2, colData=cData, design=~replicates)  #batch+


# ---- ACTIVATED for non-fake replicates ----
replicates = as.factor(c("untreated","untreated", "p601_d3","p604_d3", "p601_d8","p604_d8",   "p601_d3","p604_d3","p601_d8","p604_d8",  "p600_d3","p605_d3","p600_d8", "p605_d8",  "p600_d3","p605_d3","p600_d8","p605_d8",  "p602_d8","p602_d8",  "p600and601_d12", "p604and605_d12","p600and601_d9","p604and605_d9",     "p602_d3","p602_d3"))

batch = as.factor(c("200420", "200420", "190927", "190927",    "190927", "190927", "190228", "190228",    "190228", "190228", "191008", "191008",    "191008", "191008", "190228", "190228",     "190228", "190228", "200817", "200817",       "200420", "200420", "200817", "200817",      "210302", "210302"))

ids = as.factor(c("untreated_DonorI","untreated_DonorII", "p601_d3_DonorII","p604_d3_DonorII", "p601_d8_DonorII","p604_d8_DonorII",   "p601_d3_DonorI","p604_d3_DonorI","p601_d8_DonorI","p604_d8_DonorI",  "p600_d3_DonorII","p605_d3_DonorII","p600_d8_DonorII", "p605_d8_DonorII",  "p600_d3_DonorI","p605_d3_DonorI","p600_d8_DonorI","p605_d8_DonorI",  "p602_d8_DonorII","p602_d8_DonorI",  "p600and601_d12_DonorI", "p604and605_d12_DonorI","p600and601_d9_DonorII","p604and605_d9_DonorII",        "p602_d3_DonorI","p602_d3_DonorII"))

donor = as.factor(c("DonorI","DonorII", "DonorII","DonorII", "DonorII","DonorII",   "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorII","DonorII", "DonorII",  "DonorI","DonorI","DonorI","DonorI",  "DonorII","DonorI",  "DonorI", "DonorI","DonorII","DonorII",        "DonorI","DonorII"))

cData = data.frame(row.names=colnames(reordered.raw), replicates=replicates, donor=donor, batch=batch, ids=ids)
dds<-DESeqDataSetFromMatrix(countData=reordered.raw, colData=cData, design=~batch+replicates) # ERROR due to the factor 'batch'
dds<-DESeqDataSetFromMatrix(countData=reordered.raw, colData=cData, design=~donor+replicates)
```


## 6, convert bam to bigwig using deepTools by feeding inverse of DESeq’s size Factor
You can read the details laid out by ATpoint about how to use a scale factor with deepTools. He specified it for ATAC-seq/ChIP-seq, but the principles are the same for RNA-seq: calculate a scaling factor with DESeq2 and supply the inverse (!) to bamCoverage --scaleFactor.

https://www.biostars.org/p/317023/

https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

```sh
#Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
sizeFactors(dds)
View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2150570000
    --ignoreForNormalization chrX
    --extendReads
```
Sure, you can use the DESeq2 scale factor. I don't recall whether DESeq2 is dividing by the scale factor or multiplying by it. If it's doing the latter then you don't need to invert (just compare a normalized and raw count in DESeq2 to be sure).
```sh
--skipNonCoveredRegions --binSize 10 --scaleFactor 1/DESeq's size factor
```
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3941356

BamCoverage from deepTools (version 3.0.2) was used to generate bigwig tracks with parameters “--skipNonCoveredRegions --binSize 10 --scaleFactor 1/DESeq’s size Factor”.

```sh
#> head(reordered.raw,0)
 [1] V_8_0_mock_DonorI             V_8_0_mock_DonorII           
 [3] V_8_1_5_p601_d3_DonorII       V_8_1_5_p604_d3_DonorII      
 [5] V_8_1_5_p601_d8_DonorII       V_8_1_5_p604_d8_DonorII      
 [7] V_8_1_6_p601_d3_DonorI        V_8_1_6_p604_d3_DonorI       
 [9] V_8_1_6_p601_d8_DonorI        V_8_1_6_p604_d8_DonorI       
[11] V_8_2_3_p600_d3_DonorII       V_8_2_3_p605_d3_DonorII      
[13] V_8_2_3_p600_d8_DonorII       V_8_2_3_p605_d8_DonorII      
[15] V_8_2_4_p600_d3_DonorI        V_8_2_4_p605_d3_DonorI       
[17] V_8_2_4_p600_d8_DonorI        V_8_2_4_p605_d8_DonorI       
[19] V_8_4_1_p602_d8_DonorII       V_8_4_1_p602_d8_DonorI       
[21] V_8_3_1_p600and601_d12_DonorI V_8_3_1_p604and605_d12_DonorI
[23] V_8_3_2_p600and601_d9_DonorII V_8_3_2_p604and605_d9_DonorII
[25] V_8_4_2_p602_d3_DonorI        V_8_4_2_p602_d3_DonorII  

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

> sizeFactors(dds)
NULL
> dds <- estimateSizeFactors(dds)
> sizeFactors(dds)
            V_8_0_mock_DonorI            V_8_0_mock_DonorII 
                    0.5695507                     0.7902143 
      V_8_1_5_p601_d3_DonorII       V_8_1_5_p604_d3_DonorII 
                    0.9260477                     0.9326776 
                    
      V_8_1_5_p601_d8_DonorII       V_8_1_5_p604_d8_DonorII 
                    0.8262960                     0.8322237 
       V_8_1_3_p601_d3_DonorI        V_8_1_3_p604_d3_DonorI 
                    1.2701328                     1.2327052 
                    #3
       V_8_1_3_p601_d8_DonorI        V_8_1_3_p604_d8_DonorI 
                    1.1485877                     1.2048921 
      V_8_2_3_p600_d3_DonorII       V_8_2_3_p605_d3_DonorII 
                    0.8582589                     0.8642480 
                    #4
      V_8_2_3_p600_d8_DonorII       V_8_2_3_p605_d8_DonorII 
                    0.6629930                     0.8707503 
       V_8_2_1_p600_d3_DonorI        V_8_2_1_p605_d3_DonorI 
                    1.2831708                     1.3655531 
                    #5
       V_8_2_1_p600_d8_DonorI        V_8_2_1_p605_d8_DonorI 
                    1.1263249                     1.6911329 
      V_8_4_1_p602_d8_DonorII        V_8_4_2_p602_d8_DonorI 
                    1.4155030                     1.0738455 
                    #6
V_8_3_1_p600and601_d12_DonorI V_8_3_1_p604and605_d12_DonorI 
                    0.9139917                     0.7865048 
V_8_3_2_p600and601_d9_DonorII V_8_3_2_p604and605_d9_DonorII 
                    1.2417217                     1.0156636 
-->
> sizeFactors(dds)
            V_8_0_mock_DonorI            V_8_0_mock_DonorII 
                    0.5697385                     0.7878847 
      V_8_1_5_p601_d3_DonorII       V_8_1_5_p604_d3_DonorII 
                    0.9301130                     0.9327397 
      V_8_1_5_p601_d8_DonorII       V_8_1_5_p604_d8_DonorII 
                    0.8267532                     0.8305761 
       V_8_1_6_p601_d3_DonorI        V_8_1_6_p604_d3_DonorI 
                    1.2743873                     1.2373613 
       V_8_1_6_p601_d8_DonorI        V_8_1_6_p604_d8_DonorI 
                    1.1474853                     1.2071788 
      V_8_2_3_p600_d3_DonorII       V_8_2_3_p605_d3_DonorII 
                    0.8610875                     0.8676190 
      V_8_2_3_p600_d8_DonorII       V_8_2_3_p605_d8_DonorII 
                    0.6627211                     0.8716474 
       V_8_2_4_p600_d3_DonorI        V_8_2_4_p605_d3_DonorI 
                    1.2870415                     1.3716371 
       V_8_2_4_p600_d8_DonorI        V_8_2_4_p605_d8_DonorI 
                    1.1292211                     1.6931370 
      V_8_4_1_p602_d8_DonorII        V_8_4_1_p602_d8_DonorI 
                    1.4178388                     1.0752340 
V_8_3_1_p600and601_d12_DonorI V_8_3_1_p604and605_d12_DonorI 
                    0.9149559                     0.7891984 
V_8_3_2_p600and601_d9_DonorII V_8_3_2_p604and605_d9_DonorII 
                    1.2432898                     1.0176785 
       V_8_4_2_p602_d3_DonorI       V_8_4_2_p602_d3_DonorII 
                    0.8947267                     1.0892385 
> sizeFactors(dds)
    V_8_1_6_p601_d3_Donor1              V_8_1_5_p601_d3_Donor2 
                 2,0908056 (1,2743873)     1,5683920 (should 0,9301130) 
V_8_5_1_sTplusLT_d3_Donor1              V_8_5_2_sTplusLT_d3_Donor2 
                 0,6055365 (0,364027306)   0,5366269 (0,322601271)

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#unter konsole
for sample in V_8_0_mock_DonorI            V_8_0_mock_DonorII       V_8_1_5_p601_d3_DonorII       V_8_1_5_p604_d3_DonorII    V_8_1_5_p601_d8_DonorII       V_8_1_5_p604_d8_DonorII                   V_8_1_3_p601_d3_DonorI        V_8_1_3_p604_d3_DonorI      V_8_1_3_p601_d8_DonorI        V_8_1_3_p604_d8_DonorI     V_8_2_3_p600_d3_DonorII       V_8_2_3_p605_d3_DonorII                      V_8_2_3_p600_d8_DonorII       V_8_2_3_p605_d8_DonorII     V_8_2_1_p600_d3_DonorI        V_8_2_1_p605_d3_DonorI     V_8_2_1_p600_d8_DonorI        V_8_2_1_p605_d8_DonorI                     V_8_4_1_p602_d8_DonorII        V_8_4_2_p602_d8_DonorI     V_8_3_1_p600and601_d12_DonorI V_8_3_1_p604and605_d12_DonorI    V_8_3_2_p600and601_d9_DonorII V_8_3_2_p604and605_d9_DonorII; do
  echo "bamCoverage --bam ${sample}Aligned.sortedByCoord.out.bam -o ../bigWigs/${sample}_norm.bw --binSize 10 --scaleFactor 1/0.5695507 --effectiveGenomeSize 2864785220"
done

bamCoverage --bam ${sample}Aligned.sortedByCoord.out.bam -o ../bigWigs/${sample}_norm.bw --binSize 10 --scaleFactor 1/DESeqs_size_factor --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_0_mock_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_0_mock_DonorI_norm.bw --binSize 10 --scaleFactor 1.755769942869002 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_0_mock_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_0_mock_DonorII_norm.bw --binSize 10 --scaleFactor 1.2654795034714001 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_1_5_p601_d3_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_5_p601_d3_DonorII_norm.bw --binSize 10 --scaleFactor 1.0798579813977185 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_1_5_p604_d3_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_5_p604_d3_DonorII_norm.bw --binSize 10 --scaleFactor 1.0721818557666658 --effectiveGenomeSize 2864785220

bamCoverage --bam V_8_1_5_p601_d8_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_5_p601_d8_DonorII_norm.bw --binSize 10 --scaleFactor 1.2102200664168772 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_1_5_p604_d8_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_5_p604_d8_DonorII_norm.bw --binSize 10 --scaleFactor 1.201600002499328 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_1_3_p601_d3_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_3_p601_d3_DonorI_norm.bw --binSize 10 --scaleFactor 0.7873192472472169 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_1_3_p604_d3_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_3_p604_d3_DonorI_norm.bw --binSize 10 --scaleFactor 0.8112239649836798 --effectiveGenomeSize 2864785220
#3
bamCoverage --bam V_8_1_3_p601_d8_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_3_p601_d8_DonorI_norm.bw --binSize 10 --scaleFactor 0.870634432181365 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_1_3_p604_d8_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_1_3_p604_d8_DonorI_norm.bw --binSize 10 --scaleFactor 0.8299498353421024 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_2_3_p600_d3_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_3_p600_d3_DonorII_norm.bw --binSize 10 --scaleFactor 1.1651495836512735 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_2_3_p605_d3_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_3_p605_d3_DonorII_norm.bw --binSize 10 --scaleFactor 1.1570752839462746 --effectiveGenomeSize 2864785220
#4
bamCoverage --bam V_8_2_3_p600_d8_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_3_p600_d8_DonorII_norm.bw --binSize 10 --scaleFactor 1.508311550800687 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_2_3_p605_d8_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_3_p605_d8_DonorII_norm.bw --binSize 10 --scaleFactor 1.1484348612914632 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_2_1_p600_d3_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_1_p600_d3_DonorI_norm.bw --binSize 10 --scaleFactor 0.7793194795268097 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_2_1_p605_d3_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_1_p605_d3_DonorI_norm.bw --binSize 10 --scaleFactor 0.7323040019461711 --effectiveGenomeSize 2864785220
#5
bamCoverage --bam V_8_2_1_p600_d8_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_1_p600_d8_DonorI_norm.bw --binSize 10 --scaleFactor 0.8878432857162263 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_2_1_p605_d8_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_2_1_p605_d8_DonorI_norm.bw --binSize 10 --scaleFactor 0.5913195822752901 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_4_1_p602_d8_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_4_1_p602_d8_DonorII_norm.bw --binSize 10 --scaleFactor 0.7064626496729431 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_4_2_p602_d8_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_4_2_p602_d8_DonorI_norm.bw --binSize 10 --scaleFactor 0.9312326587018337 --effectiveGenomeSize 2864785220
#6
bamCoverage --bam V_8_3_1_p600and601_d12_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_3_1_p600and601_d12_DonorI_norm.bw --binSize 10 --scaleFactor 1.0941018392180148 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_3_1_p604and605_d12_DonorIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_3_1_p604and605_d12_DonorI_norm.bw --binSize 10 --scaleFactor 1.2714480572782263 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_3_2_p600and601_d9_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_3_2_p600and601_d9_DonorII_norm.bw --binSize 10 --scaleFactor 0.8053334334094346 --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_3_2_p604and605_d9_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_3_2_p604and605_d9_DonorII_norm.bw --binSize 10 --scaleFactor 0.9845779645937887 --effectiveGenomeSize 2864785220
#7
bamCoverage --bam V_8_4_2_p602_d3_DonorIAligned.sortedByCoord.out.bam  -o ../bigWigs/V_8_4_2_p602_d3_DonorI_norm.bw --binSize 10 --scaleFactor 0.8947267  --effectiveGenomeSize 2864785220
bamCoverage --bam V_8_4_2_p602_d3_DonorIIAligned.sortedByCoord.out.bam -o ../bigWigs/V_8_4_2_p602_d3_DonorII_norm.bw --binSize 10 --scaleFactor 1.0892385 --effectiveGenomeSize 2864785220
```


## 7, pca and heatmap before and after removing batch effects
```sh
#--go to BREAK POINT--
# sorted by sT and truncLT

#vsd <- vst(dds)
rld <- rlogTransformation(dds)

# -- before pca --
png("pca_before_batch_correction.png", 1200, 800)
plotPCA(rld, intgroup=c("replicates"))
#plotPCA(rld, intgroup = c("replicates", "batch"))
#plotPCA(rld, intgroup = c("replicates", "ids"))
#plotPCA(rld, "batch")
dev.off()


# -- before heatmap --
## generate the pairwise comparison between samples
library(gplots) 
library("RColorBrewer")
png("heatmap_before_batch_correction.png", 1200, 800)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
#paste( rld$dex, rld$cell, sep="-" )
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(replicates,batch, sep=":"))
#rownames(mat) <- colnames(mat) <- with(colData(dds),paste(replicates,ids, sep=":"))
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(13, 13))
dev.off()


# -- remove batch effect --
#show the results which delete the batches effect
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#TODO: next time using vsd
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#how-do-i-use-vst-or-rlog-data-for-differential-testing
#It uses the design formula to calculate the within-group variability (if blind=FALSE) or the across-all-samples variability (if blind=TRUE).
#- It is possible to visualize the transformed data with batch variation removed, using the removeBatchEffect function from limma.
#- This simply removes any shifts in the log2-scale expression data that can be explained by batch.
#IMPROVE: ~batch+replicates
#https://support.bioconductor.org/p/116375/
#Have a look at the manual pages of these functions. The first sentence of that for varianceStabilizingTransformation says "This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors)." For rlog, it says "This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size."
#Do try to read the documentation and a little bit about the underlying methods, you'll find that you'll be more efficient and have much more fun with the software.
mat <- assay(vsd)
mm <- model.matrix(~replicates, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
#plotPCA(vsd)

#https://www.biostars.org/p/403053/
dds <- DESeqDataSetFromMatrix(countData=counts, colData=factors, design = ~ Batch + Covariate)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$Batch)
assay(vsd) <- mat
counts_batch_corrected <- assay(vsd)


# -- after pca --
png("pca_after_batch_correction.png", 1200, 800)
#plotPCA(rld, intgroup = c("replicates", "batch"))
#plotPCA(rld, intgroup = c("replicates", "ids"))
plotPCA(rld, intgroup=c("replicates"))
dev.off()

# -- after heatmap --
## generate the pairwise comparison between samples
png("heatmap_after_batch_correction.png", 1200, 800)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(replicates,batch, sep=":"))
#rownames(mat) <- colnames(mat) <- with(colData(dds),paste(replicates,ids, sep=":"))
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(13, 13))
dev.off()
```


## 8, select the differentially expressed genes (sT)
```sh
#p600and601_d9or12
#p600_d3
#p601_d3
#p601_d8
#p602_d8
#p604and605_d9or12
#p604_d3
#p604_d8
#p605_d3
#p605_d8
#untreated
#p600_d8

#For each sample there should be a duplicate… In general we need the following comparisons:
#P600 vs untreated / mock (donor 1+2) --> p600_d3_vs_untreated, p600_d8_vs_untreated
#P601 vs untreated / mock (donor 1+2) --> p601_d3_vs_untreated, p601_d8_vs_untreated
#P604 d3 vs p601 d3 (Donor 1+2)       --> p604_d3_vs_p601_d3
#P604 d8 vs p601 d8 (donor 1+2)       --> p604_d8_vs_p601_d8
#P605 d3 vs p600 d3 (donor 1+2)       --> p605_d3_vs_p600_d3
#P605 d8 vs p600 d8 (donor 1+2)       --> p605_d8_vs_p600_d8
#P602 d8 vs p600 d8 (donor 1+2)       --> p602_d8_vs_p600_d8
#P604+605 d9/12 vs p600+601 d9/12 (donor 1+2)  --> p604and605_d9or12_vs_p600and601_d9or12

setwd("/mnt/Seagate_Corona/Data_Denise_RNASeq/results/featureCounts/degenes_2021")
#---- * to untreated ----
dds$replicates <- relevel(dds$replicates, "untreated")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p600_d3_vs_untreated","p600_d8_vs_untreated","p601_d3_vs_untreated","p601_d8_vs_untreated")

#----
dds$replicates <- relevel(dds$replicates, "p601_d3")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p604_d3_vs_p601_d3")

#----
dds$replicates <- relevel(dds$replicates, "p601_d8")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p604_d8_vs_p601_d8")

#----
dds$replicates <- relevel(dds$replicates, "p600_d3")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p602_d3_vs_p600_d3")

#----
dds$replicates <- relevel(dds$replicates, "p600_d8")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p605_d8_vs_p600_d8")

#----
dds$replicates <- relevel(dds$replicates, "p600_d8")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p605_d8_vs_p600_d8", "p602_d8_vs_p600_d8")

#----
dds$replicates <- relevel(dds$replicates, "p600and601_d9or12")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p604and605_d9or12_vs_p600and601_d9or12")


setwd("/mnt/Seagate_Corona/Data_Denise_RNASeq/results/featureCounts/degenes_DonorII")
#---- * to untreated ----
dds$replicates <- relevel(dds$replicates, "untreated_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#--> dds <- estimateSizeFactors(dds)
#--> sizeFactors(dds)
resultsNames(dds)
clist <- c("p600_d3_DonorII_vs_untreated_DonorII","p600_d8_DonorII_vs_untreated_DonorII","p601_d3_DonorII_vs_untreated_DonorII","p601_d8_DonorII_vs_untreated_DonorII")

for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}

#----
dds$replicates <- relevel(dds$replicates, "p601_d3_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p604_d3_DonorII_vs_p601_d3_DonorII")

for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}

#----
dds$replicates <- relevel(dds$replicates, "p601_d8_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p604_d8_DonorII_vs_p601_d8_DonorII")

for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}

#----
dds$replicates <- relevel(dds$replicates, "p600_d3_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p605_d3_DonorII_vs_p600_d3_DonorII")

for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}


#----
dds$replicates <- relevel(dds$replicates, "p600_d8_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p605_d8_DonorII_vs_p600_d8_DonorII")

for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}


#----
dds$replicates <- relevel(dds$replicates, "p600_d8_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p605_d8_DonorII_vs_p600_d8_DonorII", "p602_d8_DonorII_vs_p600_d8_DonorII")

for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}


#----
dds$replicates <- relevel(dds$replicates, "p600and601_d9_DonorII")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("p604and605_d9_DonorII_vs_p600and601_d9_DonorII")



# >>> math.log(3,2)  
# 1.5849625007211563
for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  geness_res <- merge(geness, res_df)
  dim(geness_res)
  write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "all.txt", sep="-"))
  #write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="-"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="-"))
}

mv p600_d3_DonorII_vs_untreated_DonorII-all.txt p600_d3_vs_untreated-all.txt
mv p600_d3_DonorII_vs_untreated_DonorII-up.txt p600_d3_vs_untreated-up.txt
mv p600_d3_DonorII_vs_untreated_DonorII-down.txt p600_d3_vs_untreated-down.txt
mv p600_d8_DonorII_vs_untreated_DonorII-all.txt p600_d8_vs_untreated-all.txt
mv p600_d8_DonorII_vs_untreated_DonorII-up.txt p600_d8_vs_untreated-up.txt
mv p600_d8_DonorII_vs_untreated_DonorII-down.txt p600_d8_vs_untreated-down.txt
mv p601_d3_DonorII_vs_untreated_DonorII-all.txt p601_d3_vs_untreated-all.txt
mv p601_d3_DonorII_vs_untreated_DonorII-up.txt p601_d3_vs_untreated-up.txt
mv p601_d3_DonorII_vs_untreated_DonorII-down.txt p601_d3_vs_untreated-down.txt
mv p601_d8_DonorII_vs_untreated_DonorII-all.txt p601_d8_vs_untreated-all.txt
mv p601_d8_DonorII_vs_untreated_DonorII-up.txt p601_d8_vs_untreated-up.txt
mv p601_d8_DonorII_vs_untreated_DonorII-down.txt p601_d8_vs_untreated-down.txt
mv p604_d3_DonorII_vs_p601_d3_DonorII-all.txt p604_d3_vs_p601_d3-all.txt
mv p604_d3_DonorII_vs_p601_d3_DonorII-up.txt p604_d3_vs_p601_d3-up.txt 
mv p604_d3_DonorII_vs_p601_d3_DonorII-down.txt p604_d3_vs_p601_d3-down.txt
mv p604_d8_DonorII_vs_p601_d8_DonorII-all.txt p604_d8_vs_p601_d8-all.txt
mv p604_d8_DonorII_vs_p601_d8_DonorII-up.txt p604_d8_vs_p601_d8-up.txt
mv p604_d8_DonorII_vs_p601_d8_DonorII-down.txt p604_d8_vs_p601_d8-down.txt
mv p605_d3_DonorII_vs_p600_d3_DonorII-all.txt p605_d3_vs_p600_d3-all.txt
mv p605_d3_DonorII_vs_p600_d3_DonorII-up.txt p605_d3_vs_p600_d3-up.txt
mv p605_d3_DonorII_vs_p600_d3_DonorII-down.txt p605_d3_vs_p600_d3-down.txt
mv p605_d8_DonorII_vs_p600_d8_DonorII-all.txt p605_d8_vs_p600_d8-all.txt
mv p605_d8_DonorII_vs_p600_d8_DonorII-up.txt p605_d8_vs_p600_d8-up.txt
mv p605_d8_DonorII_vs_p600_d8_DonorII-down.txt p605_d8_vs_p600_d8-down.txt
mv p602_d8_DonorII_vs_p600_d8_DonorII-all.txt p602_d8_vs_p600_d8-all.txt
mv p602_d8_DonorII_vs_p600_d8_DonorII-up.txt p602_d8_vs_p600_d8-up.txt
mv p602_d8_DonorII_vs_p600_d8_DonorII-down.txt p602_d8_vs_p600_d8-down.txt
mv p604and605_d9_DonorII_vs_p600and601_d9_DonorII-all.txt _p604a605_vs_p600a601-all.txt
mv p604and605_d9_DonorII_vs_p600and601_d9_DonorII-up.txt _p604a605_vs_p600a601-up.txt
mv p604and605_d9_DonorII_vs_p600and601_d9_DonorII-down.txt _p604a605_vs_p600a601-down.txt

~/Tools/csv2xls-0.4/csv_to_xls.py \
p600_d3_vs_untreated-all.txt \
p600_d3_vs_untreated-up.txt \
p600_d3_vs_untreated-down.txt \
p600_d8_vs_untreated-all.txt \
p600_d8_vs_untreated-up.txt \
p600_d8_vs_untreated-down.txt \
p601_d3_vs_untreated-all.txt \
p601_d3_vs_untreated-up.txt \
p601_d3_vs_untreated-down.txt \
p601_d8_vs_untreated-all.txt \
p601_d8_vs_untreated-up.txt \
p601_d8_vs_untreated-down.txt \
p604_d3_vs_p601_d3-all.txt \
p604_d3_vs_p601_d3-up.txt \
p604_d3_vs_p601_d3-down.txt \
p604_d8_vs_p601_d8-all.txt \
p604_d8_vs_p601_d8-up.txt \
p604_d8_vs_p601_d8-down.txt \
p605_d3_vs_p600_d3-all.txt \
p605_d3_vs_p600_d3-up.txt \
p605_d3_vs_p600_d3-down.txt \
p605_d8_vs_p600_d8-all.txt \
p605_d8_vs_p600_d8-up.txt \
p605_d8_vs_p600_d8-down.txt \
p602_d8_vs_p600_d8-all.txt \
p602_d8_vs_p600_d8-up.txt \
p602_d8_vs_p600_d8-down.txt \
_p604a605_vs_p600a601-all.txt \
_p604a605_vs_p600a601-up.txt \
_p604a605_vs_p600a601-down.txt \
-d$',' -o degenes_DonorII.xls;


~/Tools/csv2xls-0.4/csv_to_xls.py \
p602_d3_vs_p600_d3-all.txt \
p602_d3_vs_p600_d3-up.txt \
p602_d3_vs_p600_d3-down.txt \
-d$',' -o p602_d3_vs_p600_d3.xls;
```


## 9, prepare IPA tables (sT)
```sh
# First thing: commenting the 4 lines in the for-loop above.
cut -d',' -f1-1 mock_sT_d8_vs_mock_sT_d3_background.txt > symbol1.txt
cut -d',' -f1-1 sT_d3_vs_mock_sT_d3_background.txt > symbol2.txt
cut -d',' -f1-1 sT_d8_vs_mock_sT_d8_background.txt > symbol3.txt
cut -d',' -f1-1 sT_d8_vs_sT_d3_background.txt > symbol4.txt

cut -d',' -f1-7 sT_d3_vs_mock_sT_d3_background.txt > sT_d3_vs_mock_d3.txt
cut -d',' -f2-7 sT_d8_vs_mock_sT_d8_background.txt > sT_d8_vs_mock_d8.txt
cut -d',' -f2-7 mock_sT_d8_vs_mock_sT_d3_background.txt > mock_d8_vs_mock_d3.txt
cut -d',' -f2-7 sT_d8_vs_sT_d3_background.txt > sT_d8_vs_sT_d3.txt
paste -d',' sT_d3_vs_mock_d3.txt sT_d8_vs_mock_d8.txt mock_d8_vs_mock_d3.txt sT_d8_vs_sT_d3.txt > sT_degenes_IPA.txt

sT_degenes_IPA <- read.csv("sT_degenes_IPA.txt", row.names=1) 
# modify the headers
colnames(sT_degenes_IPA)<- c("sT_d3_vs_mock_d3_baseMean","sT_d3_vs_mock_d3_log2FoldChange","sT_d3_vs_mock_d3_lfcSE","sT_d3_vs_mock_d3_stat","sT_d3_vs_mock_d3_pvalue","sT_d3_vs_mock_d3_padj","sT_d8_vs_mock_d8_baseMean","sT_d8_vs_mock_d8_log2FoldChange","sT_d8_vs_mock_d8_lfcSE","sT_d8_vs_mock_d8_stat","sT_d8_vs_mock_d8_pvalue","sT_d8_vs_mock_d8_padj","mock_d8_vs_mock_d3_baseMean","mock_d8_vs_mock_d3_log2FoldChange","mock_d8_vs_mock_d3_lfcSE","mock_d8_vs_mock_d3_stat","mock_d8_vs_mock_d3_pvalue","mock_d8_vs_mock_d3_padj","sT_d8_vs_sT_d3_baseMean","sT_d8_vs_sT_d3_log2FoldChange","sT_d8_vs_sT_d3_lfcSE","sT_d8_vs_sT_d3_stat","sT_d8_vs_sT_d3_pvalue","sT_d8_vs_sT_d3_padj")

sT_degenes_IPA$SYMBOL = rownames(sT_degenes_IPA)
identical(rownames(sT_degenes_IPA), rownames(geness))
geness_sT_degenes_IPA <- merge(geness, sT_degenes_IPA)
write.csv(as.data.frame(geness_sT_degenes_IPA[order(geness_sT_degenes_IPA$sT_d8_vs_mock_d8_padj),]), file = "geness_sT_degenes_IPA.txt")
~/Tools/csv2xls-0.4/csv_to_xls.py geness_sT_degenes_IPA.txt -d$',' -o geness_sT_degenes_IPA.xls;

#under DIR degenes under KONSOLE
for comp in mock_sT_d8_vs_mock_sT_d3 sT_d3_vs_mock_sT_d3 sT_d8_vs_mock_sT_d8 sT_d8_vs_sT_d3; do \
for comp in sT_mock_d3_vs_untreated sT_mock_d8_vs_untreated LTtr_mock_d3_vs_untreated LTtr_mock_d8_vs_untreated; do \
for comp in p602_d8_vs_p600_d8 p604and605_d9ord12_vs_p600and601_d9ord12; do \
mkdir ${comp}_output; \
cut -d',' -f2- ${comp}_up.txt > ${comp}_output/upregulated_filtered; \
cut -d',' -f2- ${comp}_down.txt > ${comp}_output/downregulated_filtered; \
cut -d',' -f2- ${comp}_background.txt > ${comp}_output/background; \
cd ${comp}_output; \
~/Tools/csv2xls-0.4/csv_to_xls.py upregulated_filtered downregulated_filtered background -d$',' -o ../${comp}_degenes.xls; \
cd ..; \
done
```


## 10, KEGG pathway enrichments (sT)
```sh
##--- load the temporary results and pathways_KEGG calculation (https://yulab-smu.github.io/clusterProfiler-book/chapter6.html) ----
# under CONSOLE
# perform the GAMOLA2-annotation with “/media/jhuang/Elements/Data_Tam_RNASeq/run_with_gamola2.sh”
mkdir pathways_KEGG

#--continue from BREAK POINT--
##
#source("https://bioconductor.org/biocLite.R") 
#biocLite("org.Hs.eg.db")
#biocLite("AnnotationDbi")
library("clusterProfiler")
library("ReactomePA")
library("org.Hs.eg.db")
setwd("~/DATA/Data_Denise_RNASeq/results/featureCounts/pathways_KEGG")

for sample in mock_sT_d8_vs_mock_sT_d3 sT_d3_vs_mock_sT_d3 sT_d8_vs_mock_sT_d8 sT_d8_vs_sT_d3; do
for sample in sT_mock_d3_vs_untreated sT_mock_d8_vs_untreated LTtr_mock_d3_vs_untreated LTtr_mock_d8_vs_untreated; do \
for sample in p602_d8_vs_p600_d8 p604and605_d9ord12_vs_p600and601_d9ord12; do \
echo "temp_up <- read.csv('../degenes/${sample}_output/upregulated_filtered', row.names=1)"
echo "temp_down <- read.csv('../degenes/${sample}_output/downregulated_filtered', row.names=1)"
echo "${sample}_sig <- rbind(temp_up, temp_down)"
echo "${sample}_KEGG <- enrichKEGG(${sample}_sig\$ENTREZID)"
echo "write.table(as.data.frame(${sample}_KEGG), file = '${sample}_KEGG.txt', sep = '\t', row.names = FALSE)"
done
# copy the generated codes above to R-environment.

png("pathways_KEGG.png",width=1200, height=800)
merged_list <- merge_result(list('p602_d8_vs_p600_d8'=p602_d8_vs_p600_d8_KEGG, 'p604and605_d9ord12_vs_p600and601_d9ord12'=p604and605_d9ord12_vs_p600and601_d9ord12_KEGG))
#merged_list <- merge_result(list('sT_mock_d3-untreated'=sT_mock_d3_vs_untreated_KEGG, 'sT_mock_d8-untreated'=sT_mock_d8_vs_untreated_KEGG, 'LTtr_mock_d3-untreated'=LTtr_mock_d3_vs_untreated_KEGG, 'LTtr_mock_d8-untreated'=LTtr_mock_d8_vs_untreated_KEGG))
dotplot(merged_list, showCategory=20)
dev.off()

# under CONSOLE
cd pathways_KEGG
~/Tools/csv2xls-0.4/csv_to_xls.py sT_d3_vs_mock_sT_d3_KEGG.txt sT_d8_vs_mock_sT_d8_KEGG.txt mock_sT_d8_vs_mock_sT_d3_KEGG.txt sT_d8_vs_sT_d3_KEGG.txt -d$'\t' -o pathways_KEGG.xls
~/Tools/csv2xls-0.4/csv_to_xls.py sT_mock_d3_vs_untreated_KEGG.txt sT_mock_d8_vs_untreated_KEGG.txt LTtr_mock_d3_vs_untreated_KEGG.txt LTtr_mock_d8_vs_untreated_KEGG.txt -d$'\t' -o pathways_KEGG.xls
~/Tools/csv2xls-0.4/csv_to_xls.py p602_d8_vs_p600_d8_KEGG.txt p604and605_d9ord12_vs_p600and601_d9ord12_KEGG.txt -d$'\t' -o pathways_KEGG.xls
```


## 11, MSigDB pathway enrichments (sT)
```sh
##--- load the temporary results and pathways_MSigDB calculation on the sT-RNASeq (https://yulab-smu.github.io/clusterProfiler-book/chapter7.html) ----
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
#https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
#mkdir pathways_MSigDB
setwd("~/DATA/Data_Denise_RNASeq/results/featureCounts/pathways_MSigDB")

#gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
H_all <- read.gmt("~/REFs/MSigDB/h.all.v7.0.entrez.gmt")
C1_all <- read.gmt("~/REFs/MSigDB/c1.all.v7.0.entrez.gmt")
C2_all <- read.gmt("~/REFs/MSigDB/c2.all.v7.0.entrez.gmt")
C3_all <- read.gmt("~/REFs/MSigDB/c3.all.v7.0.entrez.gmt")
C4_all <- read.gmt("~/REFs/MSigDB/c4.all.v7.0.entrez.gmt")
C5_all <- read.gmt("~/REFs/MSigDB/c5.all.v7.0.entrez.gmt")
C6_all <- read.gmt("~/REFs/MSigDB/c6.all.v7.0.entrez.gmt")
C7_all <- read.gmt("~/REFs/MSigDB/c7.all.v7.0.entrez.gmt")

for sample in mock_sT_d8_vs_mock_sT_d3 sT_d3_vs_mock_sT_d3 sT_d8_vs_mock_sT_d8 sT_d8_vs_sT_d3; do
for sample in p602_d8_vs_p600_d8 p604and605_d9ord12_vs_p600and601_d9ord12; do \
echo "temp_up <- read.csv('../degenes/${sample}_output/upregulated_filtered', row.names=1)"
echo "temp_down <- read.csv('../degenes/${sample}_output/downregulated_filtered', row.names=1)"
echo "${sample}_sig <- rbind(temp_up, temp_down)"
echo "${sample}_MSigDB <- enricher(${sample}_sig\$ENTREZID, TERM2GENE=C7_all)"
echo "write.table(as.data.frame(${sample}_MSigDB), file = '${sample}_MSigDB.txt', sep = '\t', row.names = FALSE)"
done
# copy the generated codes above to R-environment.

# under CONSOLE
cd pathways_MSigDB
~/Tools/csv2xls-0.4/csv_to_xls.py p602_d8_vs_p600_d8_MSigDB.txt p604and605_d9ord12_vs_p600and601_d9ord12_MSigDB.txt -d$'\t' -o pathways_MSigDB_C7_all.xls
```


## 11.1, further databases for enrichments
```sh
#PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System
http://www.pantherdb.org/help/PANTHERhelp.jsp

#Phylogenetic classification of proteins encoded in complete genomes (COG)
https://www.ncbi.nlm.nih.gov/research/cog-project/

#Rfam+Pfam
http://xfam.org/
```


## 12, clustering the genes and draw heatmap (sT)
```sh
# -- prepare all_genes --
#TABOO: assay(rld) are now after batch effect removal with limma package. It is not allowed to rewrite the data with the command. rld <- rlogTransformation(dds)
RNASeq.NoCellLine <- assay(rld)
# reorder the columns
#colnames(RNASeq.NoCellLine) = c("", "")
#col.order <-  c("mock_sT_d3", "mock_sT_d3_r2", "mock_sT_d8", "mock_sT_d8_r2", "sT_d3", "sT_d3_r2", "sT_d8", "sT_d8_r2")
col.order <- c("V_8_4_1_p602_d8_DonorII","V_8_4_2_p602_d8_DonorI","V_8_2_3_p600_d8_DonorII", "V_8_2_1_p600_d8_DonorI",   "V_8_3_2_p604and605_d9_DonorII","V_8_3_1_p604and605_d12_DonorI",  "V_8_3_2_p600and601_d9_DonorII", "V_8_3_1_p600and601_d12_DonorI")
RNASeq.NoCellLine <- RNASeq.NoCellLine[,col.order]

#Option1: filter to genes with overall lfc > 2          #612
RNASeq.NoCellLine_  <- RNASeq.NoCellLine[apply(RNASeq.NoCellLine,1,function(x){max(x)-min(x)})>=2,]

#Option2: not using the automatical comparing between max and min, rather than using manually selected from DESeq2 between conditions.
cluster_ids_before <- read.csv("../../cluster_ids.csv")
geness_before <- select(org.Hs.eg.db, keys = as.vector(cluster_ids_before$ensembl_gene_id), keytype = "ENSEMBL", columns = c("SYMBOL"))  #1129
geness_before$SYMBOL
cluster1_pathways <- enrichKEGG(geness_before$SYMBOL, organism= 'hsa')
intersected_genes <- intersect(geness_before$SYMBOL, rownames(RNASeq.NoCellLine_))   #1076
RNASeq.NoCellLine_  <- RNASeq.NoCellLine[intersected_genes,]

#Option3: as paper described, A heatmap showing expression values of all DEGs which are significant between any pair conditions.
all_genes <- c(rownames(mock_sT_d8_vs_mock_sT_d3_sig),rownames(sT_d3_vs_mock_sT_d3_sig),rownames(sT_d8_vs_mock_sT_d8_sig),rownames(sT_d8_vs_sT_d3_sig))     #873
all_genes <- unique(all_genes)   #663
#all_genes2 <- c(rownames(WAC_vs_mock_sig),rownames(WAP_vs_mock_sig),rownames(WAC_vs_WAP_sig))   #3917
#all_genes2 <- unique(all_genes2)   #2608
#intersected_genes <- intersect(all_genes, all_genes2)  # 2608
#RNASeq.NoCellLine <- read.csv(file ="gene_expression_keeping_replicates.txt", row.names=1)
RNASeq.NoCellLine_  <- RNASeq.NoCellLine[all_genes,]
write.csv(as.data.frame(RNASeq.NoCellLine_), file ="gene_expression_keeping_replicates.txt")

RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, mock_sT_d3 = rowMeans(RNASeq.NoCellLine_[, 1:2]))
RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, mock_sT_d8 = rowMeans(RNASeq.NoCellLine_[, 3:4]))
RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, sT_d3 = rowMeans(RNASeq.NoCellLine_[, 5:6]))
RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, sT_d8 = rowMeans(RNASeq.NoCellLine_[, 7:8]))
RNASeq.NoCellLine_ <- RNASeq.NoCellLine_[,c(-1:-8)]        #663x4
#RNASeq.NoCellLine__ <- read.csv(file ="gene_expression_keeping_replicates.txt", row.names=1)
write.csv(as.data.frame(RNASeq.NoCellLine_ ), file ="gene_expression_merging_replicates.txt")


# -- clustering the genes and draw heatmap --
#clustering methods: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).  pearson or spearman
datamat = RNASeq.NoCellLine_
hr <- hclust(as.dist(1-cor(t(datamat), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(datamat, method="spearman")), method="complete")
mycl = cutree(hr, h=max(hr$height)/1.05)
mycol = c("YELLOW", "DARKBLUE", "DARKORANGE", "DARKMAGENTA", "DARKCYAN", "DARKRED",  "MAROON", "DARKGREEN", "LIGHTBLUE", "PINK", "MAGENTA", "LIGHTCYAN","LIGHTGREEN", "BLUE", "ORANGE", "CYAN", "RED", "GREEN");

mycol = mycol[as.vector(mycl)]
#sampleCols <- rep('GREY',ncol(RNASeq.NoCellLine_))
#names(sampleCols) <- c("mock_r1", "mock_r2", "mock_r3", "mock_r4", "WAP_r1", "WAP_r2",  "WAP_r3", "WAP_r4", "WAC_r1","WAC_r2")
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,4)=='mock'] <- 'GREY'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='WAP'] <- 'RED'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='dM_'] <- 'CYAN'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='dP_'] <- 'BLUE'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='WAC'] <- 'GREEN'
png("DEGs_heatmap.png", width=900, height=1000)
heatmap.2(as.matrix(datamat),Rowv=as.dendrogram(hr),Colv = NA, dendrogram = 'row',
            scale='row',trace='none',col=bluered(75), 
            RowSideColors = mycol, labRow="", srtCol=15)
#heatmap.2(datamat,Rowv=as.dendrogram(hr),Colv=as.dendrogram(hc),
            scale='row',trace='none',col=bluered(75),
            RowSideColors = mycol, labRow="",
            ColSideColors = sampleCols,
            main=sprintf('%s', "RNA-Seq for all DEGs"), srtCol=30)
heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none", margin=c(5,5), sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none', labRow="", srtCol=30, lhei=c(0.1,2))
#heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none",, sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none', cexRow=1.4, lhei=c(0.25,5))
#heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none", margin=c(10,10), sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none')
#heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none", margin=c(10,10), sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none', cexRow=1.2, lhei=c(0.1,5))
dev.off()


#-- cluster members --
write.csv(names(subset(mycl, mycl == '1')),file='cluster1_YELLOW.txt')
write.csv(names(subset(mycl, mycl == '2')),file='cluster2_DARKBLUE.txt') 
write.csv(names(subset(mycl, mycl == '3')),file='cluster3_DARKORANGE.txt')  
write.csv(names(subset(mycl, mycl == '4')),file='cluster4_DARKMAGENTA.txt') 
#~/Tools/csv2xls-0.4/csv_to_xls.py cluster*.txt -d',' -o genelist_clusters.xls


#-- pathway plot --
library("clusterProfiler")
library("ReactomePA")
#The cutoff of pathway enrichment is padj <=  0.001.
for sample_id in cluster1_YELLOW cluster2_DARKBLUE cluster3_DARKORANGE cluster4_DARKMAGENTA; do \
echo "${sample_id}_genes <- read.csv(\"${sample_id}.txt\")"; \
echo "geness <- select(org.Hs.eg.db, keys = as.vector(${sample_id}_genes\$x), keytype = \"SYMBOL\", columns = c(\"ENTREZID\"))"; \
echo "${sample_id}_pathways <- enrichKEGG(geness\$ENTREZID, organism= 'hsa', pvalueCutoff =  0.001, pAdjustMethod=\"BH\")"; \
echo "write.table(as.data.frame(${sample_id}_pathways), file = \"pathway_${sample_id}.txt\", sep = \"\t\", row.names = FALSE)"; \
done
png(file= 'pathways_cluster.png', width=900, height=1040)
merged_list <- merge_result(list(YELLOW=cluster1_YELLOW_pathways, DARKBLUE=cluster2_DARKBLUE_pathways, DARKORANGE=cluster3_DARKORANGE_pathways, DARKMAGENTA=cluster4_DARKMAGENTA_pathways)) 
#merge_result(list(YELLOW=cluster1_pathways, BLUE=cluster2_pathways, ORANGE=cluster3_pathways)) %>% plot(., showCategory=30, font.size = 16)
dotplot(merged_list, showCategory=30, font.size = 16)
#dotplot(merged_list, showCategory=30,)
dev.off()
#~/Tools/csv2xls-0.4/csv_to_xls.py pathway_cluster*.txt -d$'\t' -o pathways_clusters.xls


# under CONSOLE
mv gene_expression_* ../degenes/
mv DEGs_heatmap.png ../degenes/
mv cluster*.txt ../degenes/
mv genelist_clusters.xls ../degenes/
```


## 13(optional), repeat steps 8-12 for truncLT
```sh
##### STEP2: select the differentially expressed genes #####

#mock_truncLT_d3 mock_truncLT_d8 truncLT_d3 truncLT_d8
setwd("~/DATA/Data_Denise_RNASeq/results/featureCounts/degenes")
#---- relevel to mock ----
dds$replicates <- relevel(dds$replicates, "mock_truncLT_d3")
dds = DESeq(dds, betaPrior=FALSE)
resultsNames(dds)
clist <- c("mock_truncLT_d8_vs_mock_truncLT_d3", "truncLT_d3_vs_mock_truncLT_d3")

dds$replicates <- relevel(dds$replicates, "mock_truncLT_d8")
dds = DESeq(dds, betaPrior=FALSE)
clist <- c("truncLT_d8_vs_mock_truncLT_d8")

dds$replicates <- relevel(dds$replicates, "truncLT_d3")
dds = DESeq(dds, betaPrior=FALSE)
clist <- c("truncLT_d8_vs_truncLT_d3")

# >>> math.log(3,2) 
# 1.5849625007211563
for (i in clist) {
  contrast = paste("replicates", i, sep="_")
  res = results(dds, name=contrast)
  #res <- res[!is.na(res$log2FoldChange),]
  geness <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(res), keytype = "SYMBOL", columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
  geness <- geness[!duplicated(geness$SYMBOL), ]
  res$SYMBOL = rownames(res)
  rownames(geness) <- geness$SYMBOL
  identical(rownames(res), rownames(geness))
  res_df <- as.data.frame(res)
  #geness_res <- merge(geness, res_df)
  #dim(geness_res)
  #write.csv(as.data.frame(geness_res[order(geness_res$pvalue),]), file = paste(i, "background.txt", sep="_"))
  write.csv(res_df, file = paste(i, "background.txt", sep="_"))
  up <- subset(geness_res, padj<=0.05 & log2FoldChange>=2)
  down <- subset(geness_res, padj<=0.05 & log2FoldChange<=-2)
  write.csv(as.data.frame(up[order(up$log2FoldChange,decreasing=TRUE),]), file = paste(i, "up.txt", sep="_"))
  write.csv(as.data.frame(down[order(abs(down$log2FoldChange),decreasing=TRUE),]), file = paste(i, "down.txt", sep="_"))
}

############### prepare IPA tables ###############
# First thing: commenting the 4 lines in the for-loop above.
cut -d',' -f1-1 mock_truncLT_d8_vs_mock_truncLT_d3_background.txt > symbol1.txt
cut -d',' -f1-1 truncLT_d3_vs_mock_truncLT_d3_background.txt > symbol2.txt
cut -d',' -f1-1 truncLT_d8_vs_mock_truncLT_d8_background.txt > symbol3.txt
cut -d',' -f1-1 truncLT_d8_vs_truncLT_d3_background.txt > symbol4.txt

cut -d',' -f1-7 truncLT_d3_vs_mock_truncLT_d3_background.txt > truncLT_d3_vs_mock_d3.txt
cut -d',' -f2-7 truncLT_d8_vs_mock_truncLT_d8_background.txt > truncLT_d8_vs_mock_d8.txt
cut -d',' -f2-7 mock_truncLT_d8_vs_mock_truncLT_d3_background.txt > mock_d8_vs_mock_d3.txt
cut -d',' -f2-7 truncLT_d8_vs_truncLT_d3_background.txt > truncLT_d8_vs_truncLT_d3.txt
paste -d',' truncLT_d3_vs_mock_d3.txt truncLT_d8_vs_mock_d8.txt mock_d8_vs_mock_d3.txt truncLT_d8_vs_truncLT_d3.txt > truncLT_degenes_IPA.txt

truncLT_degenes_IPA <- read.csv("truncLT_degenes_IPA.txt", row.names=1) 
# modify the headers
colnames(truncLT_degenes_IPA)<- c("truncLT_d3_vs_mock_d3_baseMean","truncLT_d3_vs_mock_d3_log2FoldChange","truncLT_d3_vs_mock_d3_lfcSE","truncLT_d3_vs_mock_d3_stat","truncLT_d3_vs_mock_d3_pvalue","truncLT_d3_vs_mock_d3_padj","truncLT_d8_vs_mock_d8_baseMean","truncLT_d8_vs_mock_d8_log2FoldChange","truncLT_d8_vs_mock_d8_lfcSE","truncLT_d8_vs_mock_d8_stat","truncLT_d8_vs_mock_d8_pvalue","truncLT_d8_vs_mock_d8_padj","mock_d8_vs_mock_d3_baseMean","mock_d8_vs_mock_d3_log2FoldChange","mock_d8_vs_mock_d3_lfcSE","mock_d8_vs_mock_d3_stat","mock_d8_vs_mock_d3_pvalue","mock_d8_vs_mock_d3_padj","truncLT_d8_vs_truncLT_d3_baseMean","truncLT_d8_vs_truncLT_d3_log2FoldChange","truncLT_d8_vs_truncLT_d3_lfcSE","truncLT_d8_vs_truncLT_d3_stat","truncLT_d8_vs_truncLT_d3_pvalue","truncLT_d8_vs_truncLT_d3_padj")

truncLT_degenes_IPA$SYMBOL = rownames(truncLT_degenes_IPA)
identical(rownames(truncLT_degenes_IPA), rownames(geness))
geness_truncLT_degenes_IPA <- merge(geness, truncLT_degenes_IPA)
write.csv(as.data.frame(geness_truncLT_degenes_IPA[order(geness_truncLT_degenes_IPA$truncLT_d8_vs_mock_d8_padj),]), file = "geness_truncLT_degenes_IPA.txt")
~/Tools/csv2xls-0.4/csv_to_xls.py geness_truncLT_degenes_IPA.txt -d$',' -o geness_truncLT_degenes_IPA.xls;
############### prepare IPA tables END ###############



#under DIR degenes under KONSOLE
for comp in mock_truncLT_d8_vs_mock_truncLT_d3 truncLT_d3_vs_mock_truncLT_d3 truncLT_d8_vs_mock_truncLT_d8 truncLT_d8_vs_truncLT_d3; do \
mkdir ${comp}_output; \
cut -d',' -f2- ${comp}_up.txt > ${comp}_output/upregulated_filtered; \
cut -d',' -f2- ${comp}_down.txt > ${comp}_output/downregulated_filtered; \
cut -d',' -f2- ${comp}_background.txt > ${comp}_output/background; \
cd ${comp}_output; \
~/Tools/csv2xls-0.4/csv_to_xls.py upregulated_filtered downregulated_filtered background -d$',' -o ../${comp}_degenes.xls; \
cd ..; \
done
# 

##--- load the temporary results and save the pathways ----
# under CONSOLE
# perform the GAMOLA2-annotation with “/media/jhuang/Elements/Data_Tam_RNASeq/run_with_gamola2.sh”
mkdir pathways_KEGG

#--continue from BREAK POINT--
##
#source("https://bioconductor.org/biocLite.R") 
#biocLite("org.Hs.eg.db")
#biocLite("AnnotationDbi")
library("clusterProfiler")
library("ReactomePA")
library("org.Hs.eg.db")
setwd("~/DATA/Data_Denise_RNASeq/results/featureCounts/pathways_KEGG")

#-- mock_truncLT_d8_vs_mock_truncLT_d3 --
mock_truncLT_d8_vs_mock_truncLT_d3_up <- read.csv("../degenes/mock_truncLT_d8_vs_mock_truncLT_d3_output/upregulated_filtered", row.names=1)      #37
dim(mock_truncLT_d8_vs_mock_truncLT_d3_up)
mock_truncLT_d8_vs_mock_truncLT_d3_down <- read.csv("../degenes/mock_truncLT_d8_vs_mock_truncLT_d3_output/downregulated_filtered", row.names=1)  #1
dim(mock_truncLT_d8_vs_mock_truncLT_d3_down)
mock_truncLT_d8_vs_mock_truncLT_d3_sig <- rbind(mock_truncLT_d8_vs_mock_truncLT_d3_up, mock_truncLT_d8_vs_mock_truncLT_d3_down)                            #38
dim(mock_truncLT_d8_vs_mock_truncLT_d3_sig)
mock_truncLT_d8_vs_mock_truncLT_d3_KEGG <- enrichKEGG(mock_truncLT_d8_vs_mock_truncLT_d3_sig$ENTREZID)
write.table(as.data.frame(mock_truncLT_d8_vs_mock_truncLT_d3_KEGG), file = "mock_truncLT_d8_vs_mock_truncLT_d3_KEGG.txt", sep = "\t", row.names = FALSE)

#-- truncLT_d3_vs_mock_truncLT_d3 --
truncLT_d3_vs_mock_truncLT_d3_up <- read.csv("../degenes/truncLT_d3_vs_mock_truncLT_d3_output/upregulated_filtered", row.names=1)      #251
dim(truncLT_d3_vs_mock_truncLT_d3_up)
truncLT_d3_vs_mock_truncLT_d3_down <- read.csv("../degenes/truncLT_d3_vs_mock_truncLT_d3_output/downregulated_filtered", row.names=1)  #3
dim(truncLT_d3_vs_mock_truncLT_d3_down)
truncLT_d3_vs_mock_truncLT_d3_sig <- rbind(truncLT_d3_vs_mock_truncLT_d3_up, truncLT_d3_vs_mock_truncLT_d3_down)                              #254
dim(truncLT_d3_vs_mock_truncLT_d3_sig)
truncLT_d3_vs_mock_truncLT_d3_KEGG <- enrichKEGG(truncLT_d3_vs_mock_truncLT_d3_sig$ENTREZID)
write.table(as.data.frame(truncLT_d3_vs_mock_truncLT_d3_KEGG), file = "truncLT_d3_vs_mock_truncLT_d3_KEGG.txt", sep = "\t", row.names = FALSE)

#-- truncLT_d8_vs_mock_truncLT_d8 --
truncLT_d8_vs_mock_truncLT_d8_up <- read.csv("../degenes/truncLT_d8_vs_mock_truncLT_d8_output/upregulated_filtered", row.names=1)      #401
dim(truncLT_d8_vs_mock_truncLT_d8_up)
truncLT_d8_vs_mock_truncLT_d8_down <- read.csv("../degenes/truncLT_d8_vs_mock_truncLT_d8_output/downregulated_filtered", row.names=1)  #147
dim(truncLT_d8_vs_mock_truncLT_d8_down)
truncLT_d8_vs_mock_truncLT_d8_sig <- rbind(truncLT_d8_vs_mock_truncLT_d8_up, truncLT_d8_vs_mock_truncLT_d8_down)                               #548
dim(truncLT_d8_vs_mock_truncLT_d8_sig)
truncLT_d8_vs_mock_truncLT_d8_KEGG <- enrichKEGG(truncLT_d8_vs_mock_truncLT_d8_sig$ENTREZID)
write.table(as.data.frame(truncLT_d8_vs_mock_truncLT_d8_KEGG), file = "truncLT_d8_vs_mock_truncLT_d8_KEGG.txt", sep = "\t", row.names = FALSE)

#-- truncLT_d8_vs_truncLT_d3 --
truncLT_d8_vs_truncLT_d3_up <- read.csv("../degenes/truncLT_d8_vs_truncLT_d3_output/upregulated_filtered", row.names=1)      #24
dim(truncLT_d8_vs_truncLT_d3_up)
truncLT_d8_vs_truncLT_d3_down <- read.csv("../degenes/truncLT_d8_vs_truncLT_d3_output/downregulated_filtered", row.names=1)  #10
dim(truncLT_d8_vs_truncLT_d3_down)
truncLT_d8_vs_truncLT_d3_sig <- rbind(truncLT_d8_vs_truncLT_d3_up, truncLT_d8_vs_truncLT_d3_down)                               #34
dim(truncLT_d8_vs_truncLT_d3_sig)
truncLT_d8_vs_truncLT_d3_KEGG <- enrichKEGG(truncLT_d8_vs_truncLT_d3_sig$ENTREZID)
write.table(as.data.frame(truncLT_d8_vs_truncLT_d3_KEGG), file = "truncLT_d8_vs_truncLT_d3_KEGG.txt", sep = "\t", row.names = FALSE)

png("pathways_KEGG.png",width=1200, height=800)
merged_list <- merge_result(list('truncLT_d3-mock_truncLT_d3'=truncLT_d3_vs_mock_truncLT_d3_KEGG, 'truncLT_d8-mock_truncLT_d8'=truncLT_d8_vs_mock_truncLT_d8_KEGG, 'mock_truncLT_d8-mock_truncLT_d3'=mock_truncLT_d8_vs_mock_truncLT_d3_KEGG, 'truncLT_d8-truncLT_d3'=truncLT_d8_vs_truncLT_d3_KEGG))
dotplot(merged_list, showCategory=20)
dev.off()


# under CONSOLE
cd pathways_KEGG
~/Tools/csv2xls-0.4/csv_to_xls.py truncLT_d3_vs_mock_truncLT_d3_KEGG.txt truncLT_d8_vs_mock_truncLT_d8_KEGG.txt mock_truncLT_d8_vs_mock_truncLT_d3_KEGG.txt truncLT_d8_vs_truncLT_d3_KEGG.txt -d$'\t' -o pathways_KEGG.xls




###################################################################
##### STEP3: prepare all_genes #####
#rld <- rlogTransformation(dds)
RNASeq.NoCellLine <- assay(rld)
# reorder the columns
#colnames(RNASeq.NoCellLine) = c("", "")
col.order <-  c("mock_truncLT_d3", "mock_truncLT_d3_r2", "mock_truncLT_d8", "mock_truncLT_d8_r2", "truncLT_d3", "truncLT_d3_r2", "truncLT_d8", "truncLT_d8_r2")
RNASeq.NoCellLine <- RNASeq.NoCellLine[,col.order]

#Option1: filter to genes with overall lfc > 2          #2760
RNASeq.NoCellLine_  <- RNASeq.NoCellLine[apply(RNASeq.NoCellLine,1,function(x){max(x)-min(x)})>=2,]

#Option2: not using the automatical comparing between max and min, rather than using manually selected from DESeq2 between conditions.
cluster_ids_before <- read.csv("../../cluster_ids.csv")
geness_before <- select(org.Hs.eg.db, keys = as.vector(cluster_ids_before$ensembl_gene_id), keytype = "ENSEMBL", columns = c("SYMBOL"))  #1129
geness_before$SYMBOL
cluster1_pathways <- enrichKEGG(geness_before$SYMBOL, organism= 'hsa')
intersected_genes <- intersect(geness_before$SYMBOL, rownames(RNASeq.NoCellLine_))   #1076
RNASeq.NoCellLine_  <- RNASeq.NoCellLine[intersected_genes,]

#Option3: as paper described, A heatmap showing expression values of all DEGs which are significant between any pair conditions.
all_genes <- c(rownames(mock_truncLT_d8_vs_mock_truncLT_d3_sig),rownames(truncLT_d3_vs_mock_truncLT_d3_sig),rownames(truncLT_d8_vs_mock_truncLT_d8_sig),rownames(truncLT_d8_vs_truncLT_d3_sig))     #873
all_genes <- unique(all_genes)   #663
#all_genes2 <- c(rownames(WAC_vs_mock_sig),rownames(WAP_vs_mock_sig),rownames(WAC_vs_WAP_sig))   #3917
#all_genes2 <- unique(all_genes2)   #2608
#intersected_genes <- intersect(all_genes, all_genes2)  # 2608
#RNASeq.NoCellLine <- read.csv(file ="gene_expression_keeping_replicates.txt", row.names=1)
RNASeq.NoCellLine_  <- RNASeq.NoCellLine[all_genes,]
write.csv(as.data.frame(RNASeq.NoCellLine_), file ="gene_expression_keeping_replicates.txt")

RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, mock_truncLT_d3 = rowMeans(RNASeq.NoCellLine_[, 1:2]))
RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, mock_truncLT_d8 = rowMeans(RNASeq.NoCellLine_[, 3:4]))
RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, truncLT_d3 = rowMeans(RNASeq.NoCellLine_[, 5:6]))
RNASeq.NoCellLine_ <- cbind(RNASeq.NoCellLine_, truncLT_d8 = rowMeans(RNASeq.NoCellLine_[, 7:8]))
RNASeq.NoCellLine_ <- RNASeq.NoCellLine_[,c(-1:-8)]        #663x4
#RNASeq.NoCellLine__ <- read.csv(file ="gene_expression_keeping_replicates.txt", row.names=1)
write.csv(as.data.frame(RNASeq.NoCellLine_ ), file ="gene_expression_merging_replicates.txt")




######################################################################
##### STEP4: clustering the genes and draw heatmap #####
#clustering methods: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).  pearson or spearman
datamat = RNASeq.NoCellLine_
hr <- hclust(as.dist(1-cor(t(datamat), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(datamat, method="spearman")), method="complete")
mycl = cutree(hr, h=max(hr$height)/1.05)
mycol = c("YELLOW", "DARKBLUE", "DARKORANGE", "DARKMAGENTA", "DARKCYAN", "DARKRED",  "MAROON", "DARKGREEN", "LIGHTBLUE", "PINK", "MAGENTA", "LIGHTCYAN","LIGHTGREEN", "BLUE", "ORANGE", "CYAN", "RED", "GREEN");

mycol = mycol[as.vector(mycl)]
#sampleCols <- rep('GREY',ncol(RNASeq.NoCellLine_))
#names(sampleCols) <- c("mock_r1", "mock_r2", "mock_r3", "mock_r4", "WAP_r1", "WAP_r2",  "WAP_r3", "WAP_r4", "WAC_r1","WAC_r2")
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,4)=='mock'] <- 'GREY'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='WAP'] <- 'RED'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='dM_'] <- 'CYAN'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='dP_'] <- 'BLUE'
#sampleCols[substr(colnames(RNASeq.NoCellLine_),1,3)=='WAC'] <- 'GREEN'
png("DEGs_heatmap.png", width=900, height=1000)
heatmap.2(as.matrix(datamat),Rowv=as.dendrogram(hr),Colv = NA, dendrogram = 'row',
            scale='row',trace='none',col=bluered(75), 
            RowSideColors = mycol, labRow="", srtCol=15)
#heatmap.2(datamat,Rowv=as.dendrogram(hr),Colv=as.dendrogram(hc),
            scale='row',trace='none',col=bluered(75),
            RowSideColors = mycol, labRow="",
            ColSideColors = sampleCols,
            main=sprintf('%s', "RNA-Seq for all DEGs"), srtCol=30)
heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none", margin=c(5,5), sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none', labRow="", srtCol=30, lhei=c(0.1,2))
#heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none",, sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none', cexRow=1.4, lhei=c(0.25,5))
#heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none", margin=c(10,10), sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none')
#heatmap.2(datamat, Rowv=as.dendrogram(hr), col=bluered(75), scale="row", RowSideColors=mycol, trace="none", margin=c(10,10), sepwidth=c(0,0), dendrogram = 'row', Colv = 'false', density.info='none', cexRow=1.2, lhei=c(0.1,5))
dev.off()


#### cluster members #####
write.csv(names(subset(mycl, mycl == '1')),file='cluster1_YELLOW.txt')
write.csv(names(subset(mycl, mycl == '2')),file='cluster2_DARKBLUE.txt') 
write.csv(names(subset(mycl, mycl == '3')),file='cluster3_DARKORANGE.txt')  
#~/Tools/csv2xls-0.4/csv_to_xls.py cluster*.txt -d',' -o genelist_clusters.xls



#### pathway plot ####
library("clusterProfiler")
library("ReactomePA")
#The cutoff of pathway enrichment is padj <=  0.001 or default 0.05
for sample_id in cluster1_YELLOW cluster2_DARKBLUE cluster3_DARKORANGE; do \
echo "${sample_id}_genes <- read.csv(\"${sample_id}.txt\")"; \
echo "geness <- select(org.Hs.eg.db, keys = as.vector(${sample_id}_genes\$x), keytype = \"SYMBOL\", columns = c(\"ENTREZID\"))"; \
echo "${sample_id}_pathways <- enrichKEGG(geness\$ENTREZID, organism= 'hsa', pvalueCutoff =  0.05, pAdjustMethod=\"BH\")"; \
echo "write.table(as.data.frame(${sample_id}_pathways), file = \"pathway_${sample_id}.txt\", sep = \"\t\", row.names = FALSE)"; \
done
png(file= 'pathways_cluster.png', width=900, height=1040)
merged_list <- merge_result(list(YELLOW=cluster1_YELLOW_pathways, DARKBLUE=cluster2_DARKBLUE_pathways, DARKORANGE=cluster3_DARKORANGE_pathways)) 
#merge_result(list(YELLOW=cluster1_pathways, BLUE=cluster2_pathways, ORANGE=cluster3_pathways)) %>% plot(., showCategory=30, font.size = 16)
dotplot(merged_list, showCategory=30, font.size = 16)
#dotplot(merged_list, showCategory=30,)
dev.off()
#~/Tools/csv2xls-0.4/csv_to_xls.py pathway_cluster*.txt -d$'\t' -o pathways_clusters.xls


# under CONSOLE
mv gene_expression_* ../degenes/
mv DEGs_heatmap.png ../degenes/
mv cluster*.txt ../degenes/
mv genelist_clusters.xls ../degenes/
```
