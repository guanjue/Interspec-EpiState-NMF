The parameters in the `config.info.txt` file are defined as follows:

- script_dir: Specifies the directory path where the pipeline scripts are located. Example: /Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/scripts/.

- working_dir: Defines the working directory for the pipeline analysis. Example: /Users/guanjuexiang/Documents/projects/analysis/test_cormat_NMF_FDR_pipeline/.

- hg38_gene: The gene of interest in the human genome (hg38). Example: GATA1.

- mm10_gene: The corresponding gene of interest in the mouse genome (mm10). Example: Gata1.

- hg38_gene_set: The path to the BED file containing hg38 genes. Example: /Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/hg38.gene.bed.

- mm10_gene_set: The path to the BED file containing mm10 genes. Example: /Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/mm10.gene.bed.

- hg38_gene_exp_win_u and hg38_gene_exp_win_d: Specify the size of the upstream and downstream expansion (in base pairs) surrounding the transcription start site of the hg38 gene for the analysis. Both parameters are set to 50000 base pairs each.

- mm10_gene_exp_win_u and mm10_gene_exp_win_d: Similar to the hg38 parameters, these define the upstream and downstream expansion (in base pairs) surrounding the transcription start site of the mm10 gene. Both are set to 50000.

- hg38_state_set: The path to the BED file containing hg38 epigenetic state data. Example: /Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed. (These files can be generated using the `get_IDEAS_epigenetic_state_file.sh` script in the `/Path_to_Interspec-EpiState-NMF/scripts/`)

- mm10_state_set: The path to the BED file containing mm10 epigenetic state data. Example: /Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed. (These files can be generated using the `get_IDEAS_epigenetic_state_file.sh` script in the `/Path_to_Interspec-EpiState-NMF/scripts/`)

- EpigeneticState_meansignal_mat_file: Specifies the path to the file containing the mean signal matrix for epigenetic states. Example: /Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/IDEAS.EpigeneticState.mean_signal_mat.txt.

- fdr_threshold: The threshold for the False Discovery Rate, used in statistical significance testing. Set to 0.1.

- random_background_gene_num: The number of random background genes used in the analysis. Set to 100.

- NMF_component_num: The number of components used in the Non-negative Matrix Factorization analysis. User can choose the number of NMF factor based on the BIC figure that will be generated at the first time running the pipeline on the target genes. Example: 6.
