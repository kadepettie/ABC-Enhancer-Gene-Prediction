#! bin/bash 

#python refactor_json.py --input_json /srv/scratch/kmualim/jamboree_june_2020/hg38/input_files/input_data_lookup.json --gene_annot_collapsed_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.190910.ConsensusGeneBounds.TSS500bp.bed --nStrongestPeaks 150000 --peakExtendFromSummit 250 --outDir /srv/scratch/kmualim/jamboree_june_2020/hg38/input_files --macs_genome /mnt/lab_data2/kmualim/data/send_to_Kristy/hg38.chrom.sizes.abc --macs_type BAM 

python grabDownload.py --dhs ../output/bam_GRCh38_DHS.tsv --h3k27ac ../output/bam_GRCh38_H3K27ac.tsv --atac ../output/bam_GRCh38_ATAC.tsv  --dhs_fastq ../output/fastq_GRCh38_DHS.tsv --h3k27ac_fastq ../output/fastq_GRCh38_H3K27ac.tsv --atac_fastq ../output/fastq_GRCh38_ATAC.tsv --genome_assembly GRCh38 --outdir test_output --expt_file ../output/hg38_rel_experiments_20201224.tsv

#python grabDownload.py --dhs ../output/bam_hg19_DHS.tsv --h3k27ac ../output/bam_hg19_H3K27ac.tsv --atac ../output/bam_hg19_ATAC.tsv  --dhs_fastq ../output/fastq_hg19_DHS.tsv --h3k27ac_fastq ../output/fastq_hg19_H3K27ac.tsv --atac_fastq ../output/fastq_hg19_ATAC.tsv --genome_assembly hg19 --outdir test_output --expt_file ../output/hg19_rel_experiments_20201224.tsv 
