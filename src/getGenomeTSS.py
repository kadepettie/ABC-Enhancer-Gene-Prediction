#! bin/python3

import pandas as pd
import numpy as np
import argparse
import os
from tools import write_params, run_command
from neighborhoods import count_features_for_bed, count_single_feature_for_bed 


def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    epilog = ("")
    parser = argparse.ArgumentParser(description='Selects Top 2 Most Active TSS',
            epilog=epilog,
            formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--tss_file', required=required_args, help="tss isoform file")
    parser.add_argument('--dhs', required=required_args, help="Accessibility bam file")
    parser.add_argument('--h3k27ac', required=required_args, help="H3K27ac-seq bam file")
    parser.add_argument('--chrom_sizes', required=required_args, help="File listing chromosome size annotaions")
    parser.add_argument('--gene_outf', required=required_args, help="GeneList output")
    parser.add_argument('--genetss_outf', required=required_args, help="GeneListTSS output")
    parser.add_argument('--outDir', required=required_args)
    args = parser.parse_args()
    return args 

def read_tss_file(tss_file):
    """
    Reads in TSS File
    """
    tss_df = pd.read_csv(args.tss_file, sep="\t", names=['chr', 'start', 'end', 'TargetGeneTSS', 'score', 'strand', 'start_Gene', 'end_Gene', 'TargetGene','Collapsed_GeneIDs'])
    return tss_df

def create_dataframes(data, len_):
    df = pd.DataFrame()
    df['chr'] = data.iloc[:len_, 0].values
    df['start'] = data.iloc[len_:(len_*2), 0].values
    df['end'] = data.iloc[(len_*2): (len_*3), 0].values
    df['TargetGene'] = data.iloc[(len_*3): (len_*4), 0].values
    df['score'] = data.iloc[(len_*4): (len_*5), 0].values
    df['strand'] = data.iloc[(len_*5):, 0].values
    return df

def filter_promoters_by_distance(promoters):
    """
    Takes in Promoter Isoforms and returns Top 2 Promoter Isoforms based on RPM 
    IF promoter start sites are 500bp from each other, otherwise, return top promoter
    """
    # ensure that promoters are at least 500bp apart
    top_promoter = promoters.iloc[[0]]
    # if no promoter isoform exists within 500bp, just pick top promoter 
    # for now, use the distance between the promoter TSS 
    start = top_promoter['start'].values[0]
    if len(promoters.iloc[1:, :]) > 0:
        promoters['dist'] = promoters['start'] - start
        index = promoters.loc[np.abs(promoters['dist']) >= 500].index.astype('int')
        print("Index: ", index)
        if len(index) > 0:
            second_promoter = promoters.loc[[index[0]]]
            concatenated_top_promoters = pd.concat([top_promoter, second_promoter])
            top_promoter = filter_promoters_by_activity(concatenated_top_promoters)
            return top_promoter
        else:
            return top_promoter
    else:
        return top_promoter

def filter_promoters_by_activity(promoters):
    activities = promoters['PromoterActivityQuantile'].sort_values(ascending=False)
    activity_foldchange = activities.iloc[[1]].values[0]/activities.iloc[[0]].values[0]
    if activity_foldchange > 0.5: 
        return promoters 
    else:
        return promoters.iloc[[0]] 

def filter_expressed_df(expressed_tsscounts):
    gene_tss_df = None
    for gene in expressed_tsscounts['TargetGene'].drop_duplicates():
        tss1kb_file_subset = expressed_tsscounts.loc[expressed_tsscounts['TargetGene']==gene]
        sorted_tss1kb_file_subset = tss1kb_file_subset.sort_values(by=['PromoterActivityQuantile'], ascending=False)
        # ensure that distances between promoters are at least 500bp from each other
        top_two = filter_promoters_by_distance(sorted_tss1kb_file_subset)
        if gene_tss_df is None:
            gene_tss_df = top_two
        else:
            gene_tss_df = pd.concat([gene_tss_df, top_two])
    return gene_tss_df

def filter_nonexpressed_df(nonexpressed):
    final_df = None
    for gene in nonexpressed['TargetGene'].drop_duplicates():
        matched = nonexpressed['PromoterActivityQuantile'].loc[nonexpressed['TargetGene'].str.contains(str(gene).split("-")[0])]
        max_index = np.sort(matched).index.astype('int')
        max_match = nonexpressed.iloc[max_index, :]
        if final_df is None:
            final_df = max_match
        else:
            final_df = pd.concat([final_df, max_match])
        return final_df

def process_genome_tss(args):
    """
    Takes in ENSEMBL Gene List and outputs 1/2 Gene Isoforms for each gene 
    Promoter_ID = {Gene_Name}_{PromoterChr:Start-End}
    """
    os.makedirs(os.path.join(args.outDir), exist_ok=True)
    write_params(args, os.path.join(args.outDir, "params_generateTSS.txt"))
    
    filebase = str(os.path.basename(args.tss_file)).split(".")[0]
    feature_name =  ["H3K27ac", "DHS"]
    feature_files = [args.h3k27ac, args.dhs]
    features = {key:value for key, value in zip(feature_name, feature_files)}
    tss_df = read_tss_file(args.tss_file)
    tss1kb_file = args.tss_file
    genome_sizes = args.chrom_sizes
    outdir = args.outDir

    tsscounts = count_features_for_bed(tss_df, tss1kb_file, genome_sizes, features, outdir, "Genes.TSS1kb", force=True, use_fast_count=True)
        
    chrom_sizes = args.chrom_sizes
    tss_file = args.tss_file
    
    print("Finished Sorting Gene TSS File")
    # Take top 2 promoters based on counts 
    tsscounts['PromoterActivityQuantile'] = ((0.0001+tsscounts['H3K27ac.RPKM.quantile'])*(0.0001+tsscounts['DHS.RPKM.quantile'])).rank(method='average', na_option="top", ascending=True, pct=True)
    print("Looping though all genes present to select out Top Two Promoters based on RPM")
    tsscounts.to_csv(os.path.join(args.outDir, "PromoterActivityQuantile.tsv"), sep="\t", index=False)
    # filter for expressed genes 
    # This loop only needs to run on expressed genes
    expressed_tsscounts = tsscounts.loc[tsscounts['PromoterActivityQuantile']!=0.0]
    filtered_expressed_tsscounts = filter_expressed_df(expressed_tsscounts)
    
    nonexpressed_dup = tsscounts.loc[tsscounts['PromoterActivityQuantile']==0.0]
    # filter for single promotr entries 
    nonexpressed_unique = filter_nonexpressed_df(nonexpressed_dup)
    gene_tss_df = pd.concat([filtered_expressed_tsscounts, nonexpressed_unique])
   
    gene_tss_df.to_csv(os.path.join(args.outDir, "final.Gene_TSS.txt"), sep="\t", header=False, index=False)
    print("Saving Files")
    save_files(args, gene_tss_df)

def save_files(args, gene_tss_df):
    outdir = args.outDir
    gene_tss_df['start'] = gene_tss_df['start'].astype('int')
    gene_tss_df['end'] = gene_tss_df['end'].astype('int')
    unique_tss_df = gene_tss_df.drop_duplicates()
    unique_tss_df[['chr', 'start', 'end', 'TargetGeneTSS', 'score', 'strand', 'Collapsed_GeneIDs']].to_csv(os.path.join(outdir, "ENSEMBL_GeneTSS.txt.tmp"), sep="\t", index=False, header=False)
    get_unique_gene_id(outdir, unique_tss_df) 

def get_unique_gene_id(outdir, gene_tss_df):
    gene_tss_df['TargetGene'] = [str(gene).split("_")[1] for gene in gene_tss_df['TargetGeneTSS']]
    entries = gene_tss_df[['chr', 'TargetGene']].drop_duplicates()
    for chr, gene in zip(entries['chr'], entries['TargetGene']):
        matched = gene_tss_df.loc[(gene_tss_df['TargetGene']==gene)&(gene_tss_df['chr']==chr)]
        indexlist = matched.index.astype('int')
        if len(matched) > 1:
            count=1
            for index in indexlist:
                gene_tss_df.loc[index, 'TargetGene'] = gene_tss_df.loc[index, 'TargetGene']+"_TSS_{}".format(count)
                count+=1
    gene_tss_df[['chr', 'start_Gene', 'end_Gene', 'TargetGene', 'score', 'strand', 'Collapsed_GeneIDs']].to_csv(os.path.join(outdir, "ENSEMBL_Genes.txt.tmp"), sep="\t", index=False, header=False)

def concatenate_entries(args):
    data = pd.read_csv(os.path.join(args.outDir,"ENSEMBL_Genes.txt.tmp"), sep="\t", header=None)
    data_tss = pd.read_csv(os.path.join(args.outDir, "ENSEMBL_GeneTSS.txt.tmp"), sep="\t", header=None)
    print("Concatenating Entries for final output file")
    data.to_csv(os.path.join(args.outDir, args.gene_outf), sep="\t", header=False, index=False)
    data_tss.to_csv(os.path.join(args.outDir, args.genetss_outf), sep="\t", header=False, index=False)
    
def clean_up():
    print("Removing tmp files created during run")
    os.remove(os.path.join(args.outDir,"ENSEMBL_GeneTSS.txt.tmp"))
    os.remove(os.path.join(args.outDir,"ENSEMBL_Genes.txt.tmp"))
    print("Done!")

if __name__=="__main__":
    args = parseargs()
    outdir = args.outDir
    process_genome_tss(args)
    concatenate_entries(args)
    clean_up()
