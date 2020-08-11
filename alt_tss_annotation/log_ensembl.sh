#! bin/bash 

# filter for coding regions 
# starting with refGene.txt.gz 
# This is specific to refGene where Gene and GeneTSS 
CODEDIR=/users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/src/
DHS=/srv/scratch/kmualim/ENCODE_data/$2
H3K27ac=/srv/scratch/kmualim/ENCODE_data/$3
CELLTYPE=$1
CHROM=/mnt/lab_data2/kmualim/data/send_to_Kristy/hg19.chrom.sizes
OUTDIR=/mnt/lab_data3/kmualim/genome_files/ensembl_data/$CELLTYPE/
INPUT=Homo_sapiens.GRCh37.87.gtf.gz
OUTFILE=$INPUT.geneannotation.txt

bash /mnt/lab_data2/kmualim/ConvertBetweenAccessionID_Files/GrabAnnotationFromGeneGTF.sh $INPUT
sed -i 's/^/chr/' Homo_sapiens.GRCh37.87.gtf.gz.geneannotation.txt

# Output : _GENEINPUT.tsv 
python /mnt/lab_data3/kmualim/genome_files/scripts/filter_for_chr.py $OUTFILE

python /mnt/lab_data3/kmualim/genome_files/scripts/process_annotation_file_for_input.py ${OUTFILE}_GeneINPUT.tsv  ${OUTFILE}_GeneINPUT.tsv_CollapsedGeneBodies.txt

python /mnt/lab_data3/kmualim/genome_files/scripts/getTSS.py ${OUTFILE}_GeneINPUT.tsv_CollapsedGeneBodies.txt

bedtools slop -b 250 -i ${OUTFILE}_GeneINPUT.tsv_CollapsedGeneBodies.txt.TSS.tmp -g $CHROM > ${OUTFILE}.500bp.TSS.txt
bedtools slop -b 250 -i ${OUTFILE}_GeneINPUT.tsv_CollapsedGeneBodies.txt.TSS.tmp -g $CHROM| perl -lane '{print $F[0]."\t".$F[1]."\t".$F[2]."\t".$F[0].":".$F[1]."-".$F[2]."_".$F[3]."\t".$F[5]."\t".$F[4]}' > ${OUTFILE}.500bp.txt
cat ${OUTFILE}_GeneINPUT.tsv_CollapsedGeneBodies.txt | perl -lane '{print $F[1]."\t".$F[2]."\t".$F[3]."\t".$F[6]}' > ${OUTFILE}.tmp
paste -d"\t" ${OUTFILE}.500bp.txt ${OUTFILE}.tmp | sort -k1,1 -k2,2n > ${OUTFILE}_INPUT.tsv.sorted


file=$OUTFILE
python /mnt/lab_data3/kmualim/genome_files/scripts//getGenomeTSS.py --tss_file ${file}_INPUT.tsv.sorted --dhs $DHS --h3k27ac $H3K27ac --chrom_sizes $CHROM --gene_outf ${CELLTYPE}_Gene.txt --genetss_outf ${CELLTYPE}_GeneTSS.txt --outDir $OUTDIR 
