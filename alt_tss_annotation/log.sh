#!/bin/sh

# filter for coding regions 
# starting with refGene.txt.gz 
# This is specific to refGene where Gene and GeneTSS 
CELLTYPE=$1
DHS=/oak/stanford/groups/engreitz/Users/kmualim/ENCODEData/$2
H3K27ac=/oak/stanford/groups/engreitz/Users/kmualim/ENCODEData/$3
CHROM=/oak/stanford/groups/engreitz/Users/kmualim/ENCODEData/scripts/hg19.chrom.sizes
OUTDIR=/oak/stanford/groups/akundaje/kmualim/ABC_links/ENCODE/genome_tss/refseq_data/$CELLTYPE
OUTFILE=refGene.txt.gz.annotation

samtools index $DHS 
samtools index $H3K27ac
python getGenomeTSS.py --tss_file ${OUTFILE}_INPUT.tsv.sorted --dhs $DHS --h3k27ac $H3K27ac --chrom_sizes $CHROM --gene_outf ${CELLTYPE}_Gene.txt --genetss_outf ${CELLTYPE}_GeneTSS.txt --outDir $OUTDIR
