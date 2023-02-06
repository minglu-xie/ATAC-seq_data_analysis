#!/bin/bash
#################################################
#  File Name:record_atac_analysis.sh
#  Author: Minglu.Xie
#  Mail: minglu.xie.1023@student.uu.se,mingluxie@gmail.com
#  Created Time: Fri 03 Feb 2023 10:08:52 PM CET
#################################################

mkdir 01.merge
cat  *.filterBL.bed|sort -k1,1 -k2,2n > ./01.merge/merged.sort.filterBL.bed
bedtools merge -i ./01.merge/merged.sort.filterBL.bed > ./01.merge/merged.peaks.bed

#ls *.rmdup.bam
echo "Start to generate count matrix by bedtools multicov"
name=`ls *.rmdup.bam|awk '{ORS=" "}{print $0}'`
time bedtools multicov -bams $name -bed ./01.merge/merged.peaks.bed > 01.merge/merged.multicov.bed
echo "Matrix is ready!"
## create a title for the multicov.bed
ls *.rmdup.bam|awk '{ORS="\t"}{print $0}'|sed 's/$/\n/g'|sed 's/^/chr\tstart\tend\t/g' > ./01.merge/title_for_matrix


##create a group information for differential peak analysis (DESeq2)
echo "Start the differential peak analysis"
mkdir 02.diff_peak
ls *.rmdup.bam|awk -F "." '{if(NR<=2)print $1"\tknockout"}{if(NR>2)print $1"\tcontrol"}' > ./02.diff_peak/groupinfo

cd 02.diff_peak
/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript ../01diff_peak_V2.r ../01.merge/merged.multicov.bed groupinfo
cd ..
echo "Finish differential peak analysis"

##volcano plot
echo "Start to generate Volcano plot"
mkdir 03.volcano
cd 03.volcano
/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript ../02volcano.r ../02.diff_peak/all_peak_diff_mean.xls
cd ..
echo "Finish the volcano plot"

##Pay attention to the species; human or mouse or others?
##Input files are from original four filterBL.bed files.
##If you don't copy files to this folder, the output files will be write in the same folder as the input bed files.
echo "Start annoplot for Original Bed files"
mkdir 04.annotation
cd 04.annotation
cp ../*BL.bed ./
/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript ../03annoplot_original.r -s mm -f KOR1.filterBL.bed -d KOR2.filterBL.bed -t NULLR1.filterBL.bed -i NULLR2.filterBL.bed
cd ..


#Input files are from differential peak analysis results.
##If you don't copy files to this folder, the output files will be write in the same folder as the input files.
echo "Start annoplot and GO enrichment for differential peaks"
mkdir 05.enrichGO
cd 05.enrichGO
cp ../02.diff_peak/*.xls ./
/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript ../04enrichGO_diff.r -s mm -f knockout_up.xls -d knockout_down.xls
cd ..
echo "Finishws 05.enrichGO"

##Homer motif enrichment
echo "Start to find motif for the differential peaks"
mkdir 06.motif
cd 06.motif
awk 'NR>1{print $1"\t"$2"\t"$3}' ../02.diff_peak/knockout_up.xls > knockout_up.bed
awk 'NR>1{print $1"\t"$2"\t"$3}' ../02.diff_peak/knockout_down.xls > knockout_down.bed
/disk1/minglu/software/homer/bin/findMotifsGenome.pl knockout_up.bed mm10 knockout_up &
/disk1/minglu/software/homer/bin/findMotifsGenome.pl knockout_down.bed mm10 knockout_down &
