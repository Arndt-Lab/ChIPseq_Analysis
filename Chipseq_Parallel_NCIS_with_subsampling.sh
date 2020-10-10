#!/bin/bash
# see http://core.sam.pitt.edu/node/5678
#SBATCH -N 1
#SBATCH -t 1-00:00:00   # 1day
#SBATCH -c 8 # Request 4 cpus-per-task. #8
#SBATCH --mem=128g 
#SBATCH --job-name=ChIPseq



##LOAD MODULES

module load hisat2/2.1.0
module load deeptools/3.3.0
module load parallel/2017-05-22
module load gcc/8.2.0
module load bedtools/2.29.0
module load bamtools/2.5.1
module load samtools/1.9


parallel --citation will cite

#MAKE SCRIPT STOP ON ANY ERROR
#set -ue


#set aliases
dir=/bgfs/karndt/mae_Spt6-AID_EFs_NCIS_subsample
sample_index=sample_file_SPT6AID_EFs.txt
fastq_dir=/bgfs/karndt/mae_fastq_files
sac_cer_hisat2_index=/bgfs/karndt/Indexes/HISAT2_S_cerevisiae/genome_tran
spike_in_hisat2_index=/bgfs/karndt/Indexes/HISAT2_K_lactis/Klac

#SETUP FOR PARALLEL

#print the file name of the master file to for the user
echo "Processing all samples in file: ${sample_index}"

#get lists for parallel
#get file prefix list
file_prefix=`awk '{print $1}' ${sample_index}`
#get stain list
strain=`awk '{print $2}' ${sample_index}`
#get protein list
protein=`awk '{print $3}' ${sample_index}`
#get rep listfile_prefix = awk '{print $1}' $1
rep=`awk '{print $4}' ${sample_index}`
#get list of qubit readings
qubit=`awk '{print $5}' ${sample_index}`



#SETUP OUTPUT DIRECTORY STRUCTURE

#set working directory path
cd ${dir}

#remove old directories
rm -rf spkin_sam
rm -rf spkin_bam
rm -rf spkin_alignment_stats
rm -rf sc_sam
rm -rf sc_bam
rm -rf sc_alignment_stats
rm -rf defining_norm_regions_raw
rm -rf ncis_calculations
rm -rf counts
rm -rf bigwig_ncis
rm -rf bigwig_log2fc
rm -rf metaplots
rm -rf correlation

#make some directories
mkdir spkin_sam
mkdir spkin_bam
mkdir spkin_alignment_stats
mkdir sc_sam
mkdir sc_bam
mkdir sc_alignment_stats
mkdir defining_norm_regions_raw
mkdir ncis_calculations
mkdir counts
mkdir bigwig_ncis
mkdir bigwig_log2fc
mkdir metaplots
mkdir correlation


# ############
# ALIGNEMNT #
# ############


cd ${fastq_dir}

# #map to K.l. with hisat2
parallel --header : hisat2 \
	-p 8 \
	-q \
	-k 2 \
	--no-unal \
	--no-mixed \
	--no-discordant \
	--no-spliced-alignment \
	--no-softclip \
	-t \
	-x ${spike_in_hisat2_index} \
	-1 {SAMPLE}_R1_001.fastq.gz \
	-2 {SAMPLE}_R2_001.fastq.gz \
	--un-conc-gz {SAMPLE}_no_aln_to_spkin_R%.fastq.gz \
	--summary-file ${dir}/spkin_alignment_stats/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_all_stats_spkin.txt \
	--met-file ${dir}/spkin_alignment_stats/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_all_metrics_spkin.txt \
	-S ${dir}/spkin_sam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#convert sam to bam
parallel --header : samtools view -bh -@ 8 \
	-o ${dir}/spkin_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.bam \
	${dir}/spkin_sam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#sort bam
parallel --header : samtools sort -l 9 -O BAM -@ 8 -T ${dir}/spkin_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sorted \
	-o ${dir}/spkin_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sorted.bam \
	${dir}/spkin_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#index bam
parallel --header : samtools index -@ 8 ${dir}/spkin_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sorted.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}


map to S.c.
parallel --header : hisat2 \
	-p 8 \
	-q \
	-k 2 \
	--no-unal \
	--no-mixed \
	--no-discordant \
	--no-spliced-alignment \
	--no-softclip \
	-t \
	-x ${sac_cer_hisat2_index} \
	-1 {SAMPLE}_no_aln_to_spkin_R1.fastq.gz \
	-2 {SAMPLE}_no_aln_to_spkin_R2.fastq.gz \
	--un-conc-gz {SAMPLE}_no_aln_to_sc_R%.fastq.gz \
	--summary-file ${dir}/sc_alignment_stats/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_all_stats_sc.txt \
	--met-file ${dir}/sc_alignment_stats/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_all_metrics_sc.txt \
	-S ${dir}/sc_sam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#convert sam to bam	
parallel --header : samtools view -bh -@ 8 \
	-o ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.bam \
	${dir}/sc_sam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#https://www.biostars.org/p/365882/

#sort bam by name
The first sort can be omitted if the file is already name ordered
parallel --header : samtools sort -l 9 -n -O BAM -@ 8 -T ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sorted \
	-o ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sorted.bam \
	${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#remove unnecessary files and directories
rm -rf ${dir}/sc_sam

#run samtools fixmate
Add ms and MC tags for markdup to use later
parallel --header : samtools fixmate -rcm -@ 8 ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}.sorted.bam \
	${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_sorted_fixmate.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#remove unnecessary files and directories
rm ${dir}/sc_bam/*1.bam
rm ${dir}/sc_bam/*2.bam
rm ${dir}/sc_bam/*.sorted.bam


#sort bam
# Markdup needs position order
parallel --header : samtools sort -l 9 -O BAM -@ 8 -T ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_sorted_fixmate \
	-o ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_resorted_fixmate.bam \
	${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_sorted_fixmate.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#remove unnecessary files and directories
rm ${dir}/sc_bam/*_sorted_fixmate.bam


#remove duplicates
# Finally run mark duplicates to remove duplicates
parallel --header : samtools markdup -r -@ 8 ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_resorted_fixmate.bam \
	${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_sorted_no_duplicates.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

#remove unnecessary files and directories
rm ${dir}/sc_bam/*_resorted_fixmate.bam


#subsample BAM files for 1,000,000 read pairs (so 2,000,000 reads)
#code for this was found and nicely explained here: https://davemcg.github.io/post/easy-bam-downsampling/
#change "frac=2000000/total" to "frac=10000000/total" for example if you want 10 million reads (so 5 million pairs)

cd ${dir}/sc_bam

for BAM in `ls *_sorted_no_duplicates.bam | awk -F "." '{print $1}'`
do

	frac=$( samtools idxstats -@ 8 ${BAM}.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=2000000/total; if (frac > 1) {print 1} else {print frac}}' )

	samtools view -bs $frac -@ 8 ${BAM}.bam > ${BAM}_subsample.bam

done



cd ${dir}

#sort bam
parallel --header : samtools sort -l 9 -O BAM -@ 8 -T ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_sorted_no_duplicates_subsample \
	-o ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_subsampled_sorted.bam \
	${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_sorted_no_duplicates_subsample.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

rm ${dir}/sc_bam/*_sorted_no_duplicates_subsample.bam


#index bam
parallel --header : samtools index -@ 8 ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_subsampled_sorted.bam \
	::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}


########
# NCIS #
########

# percent = the percent of bins that will be used for the normalization
# so percent = 2 means that 2% of all genomic bins are used
# this means that effectively 2% of the genome is used for background normalization

percent=2


cd ${dir}

parallel --header : "bedtools coverage -counts -a /bgfs/karndt/Annotations/BED_files/S_c_genomic_bins_100_bp.bed -b ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_subsampled_sorted.bam > ${dir}/ncis_calculations/{STRAIN}_{PROTEIN}_{REP}_counts.bed" ::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}

parallel --header : "join -j2 -o1.1,1.2,1.3,1.4,2.4 ${dir}/ncis_calculations/{STRAIN}_Input_{REP}_counts.bed ${dir}/ncis_calculations/{STRAIN}_{PROTEIN}_{REP}_counts.bed > ${dir}/ncis_calculations/{STRAIN}_{PROTEIN}_{REP}_counts_plus_input.bed"  ::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep}



cd ${dir}/ncis_calculations

for file in `ls *_counts_plus_input.bed`
do

file_name=`echo ${file} | awk -F "." '{print $1}'`

percent_of_genome_covered_by_norm_regions=`wc -l ${file} \
| awk -v p=$percent '{print (($1/100)*p)}'`

echo percent of genome selected by user = $percent
echo file_length = `wc -l ${file}`
echo percent_of_genome_covered_by_norm_regions bin count = $percent_of_genome_covered_by_norm_regions

grep -v '#' ${file} \
| grep -v 'Mito' \
| awk ' $4>0.0 && $5>0.0 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" ($5/$4)}' \
> ${file_name}_norm_file.txt

r=1

for t in `seq 1 1 500`
do

old_r=$r

r=`awk -v t=$t ' $5<=t {ChIP+=$5; Input+=$4} END {print ChIP/Input}' ${file_name}_norm_file.txt`

n=`awk -v t=$t ' $5<=t {print}' ${file_name}_norm_file.txt | wc -l | awk '{print $1}'`

echo $t reads

if (( $(echo "$r >= $old_r" | bc -l) && $(echo "$n >= $percent_of_genome_covered_by_norm_regions" | bc -l) ))
then
	echo r=$r, n=$n, t=$t
	break
else
	echo r=$r, n=$n, t=$t
	continue
fi


done

echo norm_factor_r,num_bins_n,count_cuttoff_t > ${dir}/ncis_calculations/${file_name}_norm_stats.csv
echo $r,$n,$t >> ${dir}/ncis_calculations/${file_name}_norm_stats.csv

done 





cd ${dir}

#start file
echo norm_factor_r,num_bins_n,count_cuttoff_t > ${dir}/counts/all_norm_stats.csv
#compile stats
parallel --header : --keep-order "head -n2 ${dir}/ncis_calculations/{STRAIN}_{PROTEIN}_{REP}_counts_plus_input_norm_stats.csv | tail -n1" ::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep} >> ${dir}/counts/all_norm_stats.csv
 

#start file
echo NCIS > ${dir}/counts/NCISs.txt
##Get r value for all
awk -F "," ' NR>1 {print $1}' ${dir}/counts/all_norm_stats.csv >> ${dir}/counts/NCISs.txt

##Adjust r value to 1/r 
##because I expect background read levels to be inversely proportional to true signal
echo NCIS > ${dir}/counts/inverse_NCIS_head.txt
grep -v "NCIS" ${dir}/counts/NCISs.txt | awk -F "," '{print 1/$1}' > ${dir}/counts/inverse_NCIS.txt
cat ${dir}/counts/inverse_NCIS_head.txt ${dir}/counts/inverse_NCIS.txt > ${dir}/counts/inverse_NCISs.txt
rm ${dir}/counts/inverse_NCIS_head.txt ${dir}/counts/inverse_NCIS.txt

#append an additional column to the sample file
paste -d '\t' ${dir}/${sample_index} ${dir}/counts/inverse_NCISs.txt > ${dir}/counts/ncis_per_sample.txt


#makes region normalized bigwig
parallel --header : bamCoverage --bam ${dir}/sc_bam/{SAMPLE}_{STRAIN}_{PROTEIN}_{REP}_subsampled_sorted.bam \
--outFileName ${dir}/bigwig_ncis/{STRAIN}_{PROTEIN}_{REP}_ncis_normalized.bw \
--outFileFormat bigwig \
--scaleFactor {NCIS} \
--binSize 1 \
--ignoreDuplicates \
--extendReads \
--numberOfProcessors "max" \
::: ${file_prefix} :::+ ${strain} :::+ ${protein} :::+ ${rep} ::::+ ${dir}/counts/inverse_NCISs.txt





#####################
# SUMMARIZE BIGWIGS #
#####################

#mRNA genes
multiBigwigSummary BED-file \
	--BED /bgfs/karndt/Annotations/BED_files/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
	--bwfiles ${dir}/bigwig_ncis/*_ncis_normalized.bw \
	--blackListFileName /bgfs/karndt/Annotations/BED_files/blacklist.bed \
	--numberOfProcessors "max" \
	--outRawCounts ${dir}/correlation/multiBigwigSummary_Genes_BED.txt \
	--outFileName ${dir}/correlation/multiBigwigSummary_Genes_BED.npz 
	
#Normalization regions
multiBigwigSummary bins \
	--binSize 500 \
	--bwfiles ${dir}/bigwig_ncis/*_ncis_normalized.bw \
	--blackListFileName /bgfs/karndt/Annotations/BED_files/blacklist.bed \
	--numberOfProcessors "max" \
	--outRawCounts ${dir}/correlation/multiBigwigSummary_Bins_BED.txt \
	--outFileName ${dir}/correlation/multiBigwigSummary_Bins.npz


#####################
# PLOT CORRELATIONS #
#####################

# #mRNA genes
plotCorrelation \
--corData ${dir}/correlation/multiBigwigSummary_Genes_BED.npz \
--corMethod 'pearson' \
--whatToPlot 'heatmap' \
--plotNumbers \
--skipZeros \
--plotHeight 10 \
--plotWidth 10 \
--outFileCorMatrix ${dir}/correlation/multiBigwigSummary_Genes_BED_pearson_correlation_heatmap.tsv \
--plotFile ${dir}/correlation/multiBigwigSummary_Genes_BED_pearson_correlation_heatmap.png

# #Normalization regions
plotCorrelation \
--corData ${dir}/correlation/multiBigwigSummary_Bins.npz \
--corMethod 'pearson' \
--whatToPlot 'heatmap' \
--plotNumbers \
--skipZeros \
--plotHeight 10 \
--plotWidth 10 \
--outFileCorMatrix ${dir}/correlation/multiBigwigSummary_Bins_pearson_correlation_heatmap.tsv \
--plotFile ${dir}/correlation/multiBigwigSummary_Bins_pearson_correlation_heatmap.png






# METAPLOTS WHOLE GENE #


parallel --header : computeMatrix scale-regions \
	-b 500 -a 500 \
	-S ${dir}/bigwig_ncis/{STRAIN}_{PROTEIN}_1_ncis_normalized.bw \
	${dir}/bigwig_ncis/{STRAIN}_{PROTEIN}_2_ncis_normalized.bw \
	-R /bgfs/karndt/Annotations/BED_files/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
	--binSize 25 \
	--missingDataAsZero \
	--averageTypeBins mean \
	--numberOfProcessors "max" \
	-out ${dir}/metaplots/{STRAIN}_{PROTEIN}_matrix.gz \
	::: ${strain} :::+ ${protein}


parallel --header : plotProfile -m ${dir}/metaplots/{STRAIN}_{PROTEIN}_matrix.gz \
	-out ${dir}/metaplots/{STRAIN}_{PROTEIN}.png \
	--dpi 300 \
	--plotType=lines \
	--perGroup \
	--colors black grey \
	--samplesLabel "1" "2" \
	--legendLocation "upper-left" \
	--averageType "mean" \
	--plotHeight 8 \
	--plotWidth 8 \
	--plotFileFormat "png" \
	::: ${strain} :::+ ${protein}
# 


# METAPLOTS WHOLE GENE #

parallel --header : computeMatrix scale-regions \
	-b 500 -a 500 \
	-S ${dir}/bigwig_ncis/Spt6_AID_DMSO_{PROTEIN}_1_ncis_normalized.bw \
	${dir}/bigwig_ncis/Spt6_AID_DMSO_{PROTEIN}_2_ncis_normalized.bw \
	 ${dir}/bigwig_ncis/Spt6_AID_IAA_{PROTEIN}_1_ncis_normalized.bw \
	${dir}/bigwig_ncis/Spt6_AID_IAA_{PROTEIN}_2_ncis_normalized.bw \
	${dir}/bigwig_ncis/Spt6_AID_DMSO_IgG_1_ncis_normalized.bw \
	${dir}/bigwig_ncis/Spt6_AID_DMSO_IgG_2_ncis_normalized.bw \
	 ${dir}/bigwig_ncis/Spt6_AID_IAA_IgG_1_ncis_normalized.bw \
	${dir}/bigwig_ncis/Spt6_AID_IAA_IgG_2_ncis_normalized.bw \
	-R /bgfs/karndt/Annotations/BED_files/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed \
	--binSize 25 \
	--missingDataAsZero \
	--averageTypeBins mean \
	--numberOfProcessors "max" \
	-out ${dir}/metaplots/All_{PROTEIN}_matrix_all.gz \
	::: ${protein}


parallel --header : plotProfile -m ${dir}/metaplots/All_{PROTEIN}_matrix_all.gz \
	-out ${dir}/metaplots/All_{PROTEIN}_all.png \
	--dpi 300 \
	--plotType=lines \
	--perGroup \
	--colors black grey blue lightblue green lightgreen purple violet \
	--samplesLabel "POI_WT_1" "POI_WT_2" "POI_MUT_1" "POI_MUT_2" "IgG_WT_1" "IgG_WT_2" "IgG_MUT_1" "IgG_MUT_2" \
	--legendLocation "upper-left" \
	--averageType "mean" \
	--plotHeight 8 \
	--plotWidth 8 \
	--plotFileFormat "png" \
	::: ${protein}





exit