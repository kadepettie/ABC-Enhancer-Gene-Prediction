#! bin/bash 

# filter for coding regions 
# starting with refGene.txt.gz 
# This is specific to refGene where Gene and GeneTSS 
CODEDIR=/users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/
DHS=/srv/scratch/kmualim/ENCODE_data/$2
H3K27ac=/srv/scratch/kmualim/ENCODE_data/$3
CELLTYPE=$1
CHROM=/mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes
OUTDIR=/mnt/lab_data3/kmualim/genome_files/refseq_data/$CELLTYPE
file=refGene.txt.gz
OUTFILE=refGene.txt.gz.annotation

zcat ${file} | perl -lane 'if ($F[6]!=$F[1]) {print $F[2]."\t".$F[4]."\t".$F[5]."\t".$F[12]."\t".$F[11]."\t".$F[3]."\t".$F[1]}' > ${file}.coding_regions.tsv
python /mnt/lab_data3/kmualim/genome_files/scripts/filter_for_chr.py ${file}.coding_regions.tsv 

## lookup table 
zcat ${file} | perl -lane '{print $F[2]."\t".$F[4]."\t".$F[5]."\t".$F[1]}' > ${file}.lookuptable.txt
#
python /mnt/lab_data3/kmualim/genome_files/scripts/getTSS.py ${file}.coding_regions.tsv_GeneINPUT.tsv
#
### script outputs: ${file}.coding_regions_final.tsv.TSS.tmp
###
bedtools slop -i ${file}.coding_regions.tsv_GeneINPUT.tsv.TSS.tmp -b 250 -g $CHROM > refGene.tss.txt
bedtools slop -i ${file}.coding_regions.tsv_GeneINPUT.tsv.TSS.tmp -b 250 -g $CHROM | perl -lane '{print $F[0]."\t".$F[1]."\t".$F[2]."\t".$F[0].":".$F[1]."-".$F[2]."_".$F[3]."\t".$F[4]."\t".$F[5]}' > ${OUTFILE}.500bp.txt 
#
#
cat ${file}.coding_regions.tsv_GeneINPUT.tsv | perl -lane '{print $F[1]."\t".$F[2]."\t".$F[3]."\t".$F[6]}' > ${OUTFILE}.tmp
paste -d"\t" ${OUTFILE}.500bp.txt ${OUTFILE}.tmp | sort -k1,1 -k2,2n > ${OUTFILE}_INPUT.tsv.sorted

python /mnt/lab_data3/kmualim/genome_files/scripts/getGenomeTSS.py --tss_file ${OUTFILE}_INPUT.tsv.sorted --dhs $DHS --h3k27ac $H3K27ac --chrom_sizes $CHROM --gene_outf ${CELLTYPE}_Gene.txt --genetss_outf ${CELLTYPE}_GeneTSS.txt --outDir $OUTDIR
