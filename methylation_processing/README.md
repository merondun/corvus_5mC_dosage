# Working with 5mC

These scripts extract 5mC and run the analyses 

## File Descriptions

1. `Genome_Prep.sh`: Prepares the reference genome for bismark mapping.
2. `Bismark_5mC_Extraction.sh`: Extracts 5mC data from sequencing reads using Bismark.
3. `CpG_Summary.sh`: Generates summary statistics of CpG methylation.
4. `MBias_Summary.sh`: Summarizes methylation bias along read lengths in the data.
5. `Plot_QC.R`: Generates quality control plots.
6. `Overlap_With_Genes_CpGislands_Chromosomes.sh`: Determines overlaps between methylation data, genes, CpG islands, and chromosomes.
7. `Extract_5mC_Calls_Overlap_Across_Samples.R`: Extracts 5mC data and finds overlaps across samples using R.
8. `Final_Distributional_Summaries.R`: Creates final summary distributions of the data.
9. `Bootstrapped_FM.R`: Performs bootstrapped F/M methylation on the data.
10. `Bootstrapped_FM_Parallel_Submission.sh`: Submits the above script in parallel.
11. `Merge_with_chromosome_Data.sh`: Merges the methylation data with chromosome-wide data.
12. `Identify_Compensation_States.sh`: Identifies potential dosage compensation states in the data.
13. `Overlap_Gene_CpGIsland_Generate_States_Plot.R`: Overlaps gene and CpG island data and generates state plots.

