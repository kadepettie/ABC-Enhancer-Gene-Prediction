#! bin/bash 

import sys, os
import argparse 
import pandas as pd
import numpy as np

def is_unique(entries):
    a = entries.to_numpy()
    return (a[0] == a[1:]).all()

def check_x_and_y_single(entries, column):
    columns = column
    return_values = ['File accession_Accessibility']
    if not is_unique(entries[columns[0]]):
        return str(return_values[0]), str(return_values[1])
    else:
        return -1, -1

def check_x_and_y(entries, column):
    columns = column
    return_values = ['File accession_Accessibility bam', 'File accession_H3K27ac bam']
    # check for biological replicates
    if not is_unique(entries[columns[0]]):
        return str(return_values[0]), str(return_values[1])
    elif not is_unique(entries[columns[1]]):
        return str(return_values[1]), str(return_values[0])
    else:
        return -1, -1
    
# This function is used by grabUnique to treat all ambigiously paired bam files as singleend 
def findFilterFiles(outfile, outfile2):
    p = pd.read_csv(outfile, sep="\t", header=None)
    p2 = pd.read_csv(outfile2, sep="\t", header=None)

    common = p[0].loc[p[0].isin(p2[0])]
    remove_files_p2 = p2[0].loc[np.logical_not(p2[0].isin(list(common)))]
    remove_files_p2.to_csv(outfile2, sep="\t", header=False, index=False)

# This function grabs all unique entries in singleend and paireend lookup table for entry into removing duplicates function 
# It also resolves the conflict of bamfiles that are inaccurately labelled as both "pairedend" and "singleend" in ENCODE 
def grabUnique(outdir, filenames):
    for filename in filenames:
        data = pd.read_csv(os.path.join(outdir, filename), sep="\t", header=None)
        unique = data[0].drop_duplicates()
        unique.to_csv(os.path.join(outdir, "unique_"+str(filename)), sep="\t", index=False, header=False)
    # filter both files to ensure that each entry only appears once
    outfile = os.path.join(outdir, "unique_"+str(filenames[0]))
    outfile2 = os.path.join(outdir, "unique_"+str(filenames[1]))

    # find common entries
    findFilterFiles(outfile, outfile2)

def rename_bam_for_pairedend(metadata, col, paired):
    indicies = metadata['File accession_'+str(col)].loc[metadata['Run type_'+str(col)]==str(paired)].index.astype('int')
    metadata['File accession_'+str(col)+" bam"] = metadata['File accession_'+str(col)]
#    accessibility_col = metadata.columns.get_loc('File accession_'+str(col)+" bam")
    metadata.loc[indicies, 'File accession_'+str(col)+" bam"] = metadata.loc[indicies, 'File accession_'+str(col)+" bam"].astype('str') + ".nodup"
    return metadata
    
def rename_biosample_term_name(new_data):
    for biosample in new_data['Biosample term name'].drop_duplicates():
        matches = new_data.loc[new_data['Biosample term name']==biosample]
        indices = matches.index.astype('int')
        if len(matches) > 1: 
            index=1
            for i in indices[1:]:
                new_data.loc[i, 'Biosample term name'] = new_data.loc[i,'Biosample term name']+"_{}".format(index)
                index+=1
    return new_data
# This function prepares the lookup tables for input into ABC
# Where each DHS, H3K27ac bam file is appended to its corresponding celltype
def prepareLookup(args, metadata, title):
# for pairedend data accession files, change name to *.nodup.bam since duplicates need to be removed from these files
#    metadata_tmp = rename_bam_for_pairedend(metadata, 'Accessibility', "paired-ended") 
#    metadata = rename_bam_for_pairedend(metadata_tmp, 'H3K27ac', "paired-ended")
    # process biosample term name 
    metadata['Biosample term name'] = [str(i).replace(",", "").replace(" ", "_") for i in metadata['Biosample term name']]
    celltypes = metadata['Biosample term name'].drop_duplicates()
    celltypes.to_csv(os.path.join(args.outdir, "cells.txt"), sep="\t", header=False, index=False)
    metadata['File accession_Accessibility bam'] = metadata['File accession_Accessibility bam']+".bam"
    if args.h3k27ac is not None:
        metadata['File accession_H3K27ac bam'] = metadata['File accession_H3K27ac bam']+".bam"
        metadata[['Biosample term name', 'File accession_Accessibility bam', 'File accession_H3K27ac bam']].to_csv(os.path.join(args.outdir, str(title)+".tsv"), sep="\t", header=False, index=False)
        metadata = rename_biosample_term_name(metadata)
        new_data = metadata.rename(index={key:value for key, value in zip(metadata.index, metadata['Biosample term name'])})
        new_data[['File accession_Accessibility bam', 'File accession_H3K27ac bam']].to_json(os.path.join(args.outdir, str(title)+".json"), orient='index')
    else:
        metadata = rename_biosample_term_name(metadata)
        new_data = metadata.rename(index={key:value for key, value in zip(metadata.index, metadata['Biosample term name'])})
        new_data[['File accession_Accessibility bam']].to_json(os.path.join(args.outdir, str(title)+".json"), orient='index')
    return metadata

def mapExperimentToLength(dhs, dhs_fastq):
    dhs_lookup = dhs_fastq[['Experiment accession', 'Run type']].drop_duplicates()
    dhs_bam = dhs.loc[(dhs['File format']=='bam') & (dhs['Output type'].str.contains('unfiltered alignments'))]
    paired_type = []
    for name in dhs_bam['Experiment accession']:
        matched = dhs_lookup['Run type'].loc[dhs_lookup['Experiment accession']==name]
        if len(matched) >1 :
            type_ = matched.values[0]
        else:
            type_ = "single-ended"
        paired_type.append(type_)
    dhs_bam['Run type'] = paired_type
    return dhs_bam
        

# This function updates the lookup table such that entries that have bamfiles with biological replicates undergo a merging process and the collective pooled bam is now updated in the table 
def getExperimentsCombined(args,metadata, biosample_entries):
    to_combine = {}
    df = metadata
    update_lookup = {}
    for biosample, treatment, duration, modifications, categories, targets, gene_targets,bio_rep, bio_rep1 in zip(biosample_entries['Biosample term name'], biosample_entries['Biosample treatments'], biosample_entries['Biosample treatments duration'], biosample_entries['Biosample genetic modifications methods'],biosample_entries['Biosample genetic modifications categories'],biosample_entries['Biosample genetic modifications targets'], biosample_entries['Biosample genetic modifications gene targets'], biosample_entries['Biological replicate(s)_Accessibility'], biosample_entries['Biological replicate(s)_H3K27ac']):
        biosample_entry = metadata.loc[(metadata['Biosample term name']==biosample)&(metadata['Biosample treatments']==treatment)&(metadata['Biosample treatments duration']==duration)&(metadata['Biosample genetic modifications methods']==modifications)&(metadata['Biosample genetic modifications categories']==categories)&(metadata['Biosample genetic modifications targets']==targets)&(metadata['Biosample genetic modifications gene targets']==gene_targets)&(metadata['Biological replicate(s)_Accessibility']==bio_rep)&(metadata['Biological replicate(s)_H3K27ac']==bio_rep1)]
        if len(biosample_entry) > 1:
            if args.h3k27ac is not None:
                col1, col2 = check_x_and_y(biosample_entry, column=['Technical replicate(s)_Accessibility', 'Technical replicate(s)_H3K27ac'])
            else:
                col1, _ = check_x_and_y(biosample_entry, column=['Technical replicate(s)_Accessibility'])
            if col1 != -1:
                index = list(biosample_entry.index.astype('int'))
                df.loc[index[0], col1] = str(biosample_entry.loc[index[0], col1]) + "_pooled"
                df = df.drop(index[1:])
                to_combine[str(biosample).replace(",", "").replace(" ", "_")] = list(set([str(entry)+".bam" for entry in biosample_entry.loc[:,col1]])) + [str(biosample_entry.loc[index[0], col1]) + "_pooled.bam"]
            else:
                index = list(biosample_entry.index.astype('int'))
                df = df.drop(index[1:])
    return to_combine, df

def save_metadata(args, duplicates, prefix): 
    duplicates.to_csv(os.path.join(args.outdir, "{}_metadata.tsv".format(prefix)), sep="\t", index=False)
    # grab celltypes with biological replicates
    df_biological =  duplicates.loc[duplicates.duplicated(['Biosample term name'], keep=False)]
    df_biological.to_csv(os.path.join(args.outdir, "{}_Replicates_metadata.tsv".format(prefix)), sep="\t", index=False)

# This function grabs the samples that have paired ends for paired end bam processing 
# It saves pairedend bams and singleend bams for removal of duplicates 
# Treat biological replicates as separate samples 
# Pick out technical replicates for Accessibility and H3K27ac 
def obtainDuplicated(args, subset_intersected):
    # grab Biosample term name that have merged technical replicates 
    merged_technical_replicate_file = subset_intersected.loc[subset_intersected['Technical replicate(s)_Accessibility'].str.contains(',') | subset_intersected['Technical replicate(s)_H3K27ac'].str.contains(',')]


    merged_index = merged_technical_replicate_file.index.astype('int')
    total = subset_intersected.index.astype('int')
    single_index = [i for i in total if i not in merged_index]
    # grab Biosample term name that have single technical replicates 
    single_technical_replicate_file = subset_intersected.loc[single_index, :]
    # filter for Experiment accession in merged_technical_replicate_file
    to_filter_single_technical_replicate_file = single_technical_replicate_file.loc[np.logical_not(single_technical_replicate_file['Experiment accession_Accessibility'].isin(list(merged_technical_replicate_file['Experiment accession_Accessibility']))) | np.logical_not(single_technical_replicate_file['Experiment accession_H3K27ac'].isin(list(merged_technical_replicate_file['Experiment accession_H3K27ac'])))]

    # filter final merged file for these merged columns 
    duplicate_merge_columns = ['Biosample term name', 'Biosample organism', 'Biosample treatments','Biosample treatments amount', 'Biosample treatments duration','Biosample genetic modifications methods','Biosample genetic modifications categories','Biosample genetic modifications targets', 'Biosample genetic modifications gene targets', 'File assembly', 'Genome annotation', 'File format', 'File type', 'Output type', 'Lab', 'Biological replicate(s)_Accessibility', 'Biological replicate(s)_H3K27ac', 'Technical replicate(s)_Accessibility', 'Technical replicate(s)_H3K27ac']
    
    
    # Label duplicated experiment accession numbers (file entries from similar experiments) and removes duplicate entries based on Biological replicates 
    dhs_duplicates = to_filter_single_technical_replicate_file[to_filter_single_technical_replicate_file.duplicated(['Experiment accession_Accessibility'], keep=False)].drop_duplicates(duplicate_merge_columns)
    
    if args.h3k27ac is not None:
        h3k27acduplicates = to_filter_single_technical_replicate_file[to_filter_single_technical_replicate_file.duplicated(['Experiment accession_H3K27ac'], keep=False)].drop_duplicates(duplicate_merge_columns)
    
        comb_duplicates = pd.concat([dhs_duplicates, h3k27acduplicates])
        duplicates = comb_duplicates.drop_duplicates()
        duplicates["File accession_H3K27ac bam"] = [str(i).split(".bam")[0] for i in duplicates["File accession_H3K27ac"]] 
    else: 
        duplicates = dhs_duplicates.drop_duplicates()
    
    # save single duplicate file
    save_metadata(args, duplicates, "single_technical_rep")
    
    # process merged_technical_replicate file for pairedend reads
    #metadata_orig = pd.concat([merged_technical_replicate_file, duplicates])
    #save_metadata(args, metadata_orig, "single_technical_rep_merged_rep")
    
    # rename paired-end duplicates file
    paired_end_duplicates = rename_bam_for_pairedend(duplicates, 'Accessibility', "paired-ended")
    #save_metadata(args, duplicates_paired_end, "single_technical_rep_merged_rep_rename_accessibility_paired-end")
    
    paired_end_merged_technical_replicate_file = rename_bam_for_pairedend(merged_technical_replicate_file, "Accessibility", "paired-ended")
    paired_end_merged_technical_replicate_file["File accession_H3K27ac bam"] = [str(i).split(".bam")[0] for i in paired_end_merged_technical_replicate_file["File accession_H3K27ac"]]
    
    # index file to process for technical replicates that need to be merged 
    to_merge_technical_replicate_index = duplicates.index.astype('int')
    df_biological_rep = duplicates[['Biosample term name', 'Biosample treatments', 'Biosample treatments amount', 'Biosample treatments duration', 'Biosample genetic modifications methods', 'Biosample genetic modifications categories', 'Biosample genetic modifications targets', 'Biosample genetic modifications gene targets', 'Biological replicate(s)_Accessibility', 'Biological replicate(s)_H3K27ac']].drop_duplicates().dropna()

    to_combine, metadata_unique = getExperimentsCombined(args, duplicates, df_biological_rep)
    with open(os.path.join(args.outdir, "Experiments_ToCombine.txt.tmp"), "w") as f:
        for key, value in zip(to_combine.keys(), to_combine.values()):
            f.write(str(key))
            f.write("\t")
            f.write(str(list(value)))
            f.write("\n")
        f.close()
    
    metadata_orig = pd.concat([metadata_unique, paired_end_merged_technical_replicate_file])
    save_metadata(args,metadata_orig, "Finalized")
    return metadata_orig, metadata_unique

