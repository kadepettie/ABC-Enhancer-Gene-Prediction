#! bin/bash 

python refactor_json.py --input_json /srv/scratch/kmualim/jamboree_june_2020/hg38/input_files/input_data_lookup.json --gene_annot_collapsed_tss /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.190910.ConsensusGeneBounds.TSS500bp.bed --nStrongestPeaks 150000 --peakExtendFromSummit 250 --outDir /srv/scratch/kmualim/jamboree_june_2020/hg38/input_files --macs_genome /mnt/lab_data2/kmualim/data/send_to_Kristy/hg38.chrom.sizes.abc --macs_type BAM 
