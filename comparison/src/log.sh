#Rscript comparePredictionsToExperiment.R --predictions /mnt/lab_data2/kmualim/Jamboree_data/notebooks/NotebooksForComparingDifferentABCInputs/K562.EnhancerPredictions.Merged.TopPromoter.tsv --experimentalData ../data/ExperimentalData.Gasperini.FulcoNasser.191021.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt
#mv *.txt /oak/stanford/groups/akundaje/projects/ABC_links/plots/ComparingDiffTSSAnnotation/merged_top_promoter_K562_alt_tss/full/
#mv *.pdf /oak/stanford/groups/akundaje/projects/ABC_links/plots/ComparingDiffTSSAnnotation/merged_top_promoter_K562_alt_tss/full/
#
#Rscript comparePredictionsToExperiment.R --predictions /mnt/lab_data2/kmualim/Jamboree_data/notebooks/NotebooksForComparingDifferentABCInputs/K562.EnhancerPredictions.Merged.TopPromoter.tsv --experimentalData ../data/ALL_but_GASPERINI_Experimental.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt
#mv *.txt /oak/stanford/groups/akundaje/projects/ABC_links/plots/ComparingDiffTSSAnnotation/merged_top_promoter_K562_alt_tss/all_but_gasperini/
#mv *.pdf /oak/stanford/groups/akundaje/projects/ABC_links/plots/ComparingDiffTSSAnnotation/merged_top_promoter_K562_alt_tss/all_but_gasperini/

#Rscript comparePredictionsToExperiment.R --predictions /mnt/lab_data2/kmualim/distalreg_scripts/analyze_links_in_out_tads/K562_TAD_Domains/EnhancerPredictionsAllPutative_subCols_hic_kr_include.txt --experimentalData ../data/ExperimentalData.Gasperini.FulcoNasser.191021.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt

#bedtools intersect -a test.txt -b /mnt/lab_data2/kmualim/Jamboree_data/notebooks/NotebooksForComparingDifferentABCInputs/stripe_loops.txt -wa -wb > K562_ExperimentalData_EnhancerPredsAllPutative_hic_kr_include_stripe_loops.txt 

#Rscript comparePredictionsToExperiment.R --predictions --experimentalData ../data/ALL_but_GASPERINI_Experimental.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt
#Rscript comparePredictionsToExperiment.R --predictions /mnt/lab_data3/kmualim/Seventy_CellTypeFiles/alt_tss_predictions/K562_new_tss_new_code/Predictions/EnhancerPredictionsAllPutative_forComparison.txt.gz --experimentalData ../data/ALL_but_GASPERINI_Experimental.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt
#Rscript comparePredictionsToExperiment.R --predictions /mnt/lab_data2/kmualim/Jamboree_data/notebooks/NotebooksForComparingDifferentABCInputs/ThreeMethodIntersect.txt --experimentalData ../data/ExperimentalData.Gasperini.FulcoNasser.191021.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt

Rscript comparePredictionsToExperiment.R --predictions /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/comparison/example/input/pred.table.txt --experimentalData ../example/input/K562.ExperimentalData.slim.txt --plotConfig ../CRISPR/plot.config.txt --predConfig ../CRISPR/pred.config.txt --code /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/comparison/src/comparison.R --outDir test
