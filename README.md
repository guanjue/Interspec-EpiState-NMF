
```
 time bash /Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/cormat_NMF_FDR/Interspec-EpiState-NMF.sh /Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/cormat_NMF_FDR/input_files/config.info.txt 2> test_run.log.txt
```

# Interspec-EpiState-NMF      <img src="https://raw.githubusercontent.com/guanjue/public_log_descriptions/main/Interspec-EpiState-NMF/VisionBMG2_rh.png" align="right" width="120"/>


### Identifying inter species epigenetic state correlation between human and mouse genes by NMF

##
**[(1) Summary](#Summary)**<br>
#####
**[(2) Citation](#Citation)**<br>
#####
**[(3) Interspec-EpiState-NMF Overview](#Interspec-EpiState-NMF-Overview)**<br>
#####
**[(4) Requirements](#Requirements)**<br>
#####
**[(5) Installation](#Installation)**<br>
#####
**[(6) Input data](#Input-data)**<br>
#####
**[(7) Running Interspec-EpiState-NMF](#Running-Interspec-EpiState-NMF)**<br>
#####
**[(8) Output of Interspec-EpiState-NMF pipeline](#Output-of-Interspec-EpiState-NMF-pipeline)**<br>
#####
**[(9) Support](#Support)**<br>
#####
**[(10) LICENSE](#LICENSE)**<br>
#####

## Summary
Combinatorial patterns of epigenetic features reflect transcriptional states. Existing normalization approaches may distort relationships between functionally correlated features by normalizing each feature independently. We present JMnorm, a novel approach that normalizes multiple epigenetic features simultaneously by leveraging information from correlated features. We show that JMnorm-normalized data preserve cross-feature correlations and combinatorial patterns of epigenetic features across cell types, improve cross-cell type gene expression prediction models, consistency between biological replicates, and detection of epigenetic changes upon perturbations. These findings suggest that JMnorm minimizes technical noise while preserving biologically relevant relationships between features. 

## Citation
Guanjue Xiang, ..., Ross Hardison. Interspecies regulatory landscapes and elements revealed by novel joint systematic integration of human and mouse blood cell epigenomes. (2023)


## Interspec-EpiState-NMF Overview
![logo](https://raw.githubusercontent.com/guanjue/public_log_descriptions/main/Interspec-EpiState-NMF/XiangEtAl_JointHMVISION_Figures.png)
An overview of the four key steps in the JMnorm normalization procedure. (A) Step 1: orthogonal transformation. The correlated components of various epigenetic signals are transformed into mutually independent high-dimensional PCA dimensions. Each colored block on the left represents the signal vector of all epigenetic features at the nref or ntar cCRE regions in reference or target sample, respectively. Colored blocks on the right denote corresponding transformed PCA epigenetic signal matrices for reference and target samples. The yellow box in the middle represents the PCA rotation matrix learned from the reference signal matrix. (B) Step 2: cCRE clustering. Reference cCRE clusters are generated based on the reference data in PCA space with the average signal reference matrix shown as a heatmap. Target cCREs are assigned to reference clusters according to the Euclidean distances between the signal vector of the target cCRE and the average signal vectors of reference clusters in the PCA space. Within each cluster, the number of cCREs, shown as colored blocks within the insert, may vary between the reference and target samples. (C) Step 3: within-cluster normalization. Target signal matrix is normalized against the reference matrix using within-cluster quantile normalization as shown for Cluster k. (D) Step4: reconstruction of the JMnorm-normalized target signal matrix in the original signal space. The yellow box in the middle indicates the transposed PCA rotation matrix learned in the first step (panel A).

## Requirements
bedtools
R
R Packages: pheatmap

## Installation 
```
# create conda environment for JMnorm with required dependencies
conda config --add channels bioconda
conda create -n esnmf r r-pheatmap bedtools

# activate the environment
conda activate esnmf
```


## Input data
- The input files 

```
# list all input files
>>> cd /Path_to_Interspec-EpiState-NMF/input_files/
>>> ls
IDEAS.EpigeneticState.mean_signal_mat.txt
config.info.txt
S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed
S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed
hg38.gene.bed
mm10.gene.bed
```

- The `config.info.txt` includes all parameters used for the Interspec-EpiState-NMF pipeline. The details about the parameters in the `config.info.txt` file can be found [Here](https://raw.githubusercontent.com/guanjue/public_log_descriptions/main/Interspec-EpiState-NMF/parameter.details.md):
```
>>> cat /Path_to_Interspec-EpiState-NMF/input_files/config.info.txt
script_dir	/Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/scripts/
working_dir	/Users/guanjuexiang/Documents/projects/analysis/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/
#
hg38_gene	GATA1
mm10_gene	Gata1
hg38_gene_set	/Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/hg38.gene.bed
mm10_gene_set	/Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/mm10.gene.bed
hg38_gene_exp_win_u	50000
hg38_gene_exp_win_d	50000
mm10_gene_exp_win_u	50000
mm10_gene_exp_win_d	50000
#
hg38_state_set	/Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed
mm10_state_set	/Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed
EpigeneticState_meansignal_mat_file	/Users/guanjuexiang/Documents/projects/git/Interspec-EpiState-NMF/input_files/IDEAS.EpigeneticState.mean_signal_mat.txt
#
random_background_gene_num	100
fdr_threshold	0.1
NMF_component_num	6
```

- Due to the file size (~700 MB), the two epigenetic state bed files: `S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed` and `S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed` need to be generated by using the `get_IDEAS_epigenetic_state_file.sh` script. 

```
>>> # BED file containing data on Human epigenetic states across various cell types
>>> head S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed
CHR	POSst	POSed	AVE	B_B15_50	B_NC14_42	CMP_100246	ERY_S002R5	ERY_S002S3	GMP_100256	HSC_100258	MEP_Donor2596	MK_S004BTH2	MK_S00VHKH1	MONc_C0011IH1	NK_S005YG	T_CD4_S008H1	T_CD8_C0066PH1
chr1	792600	792800	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	792800	793000	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	793000	793200	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	793200	793400	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	793400	793600	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	793600	793800	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	793800	794000	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	794000	794200	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	794200	794400	0	0	2	0	0	0	0	0	0	2	0	0	2
```

```
>>> # BED file containing data on Mouse epigenetic states across various cell types
>>> head S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed
CHR	POSst	POSed	AVE	B_r1	B_r2	CMP_r1	ERY_ad_r1	ERY_ad_r2	GMP_r1	LSK_r1	MEP_r1	MK_fl_r1	MK_fl_r2	MON_r1	NK_r1	T_CD4_r1	T_CD8_r1
chr1	0	200	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	200	400	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	400	600	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	600	800	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	800	1000	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	1000	1200	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	1200	1400	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	1400	1600	0	0	0	0	0	0	0	0	0	0	0	0	0
chr1	1600	1800	0	0	0	0	0	0	0	0	0	0	0	0	0
```



## Running Interspec-EpiState-NMF
```
# activate the environment
conda activate esnmf

# mkdir and cd to working_dir
mkdir ~/projects/analysis/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/
cd ~/projects/analysis/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/

# run Interspec-EpiState-NMF pipeline
time bash /Path_to_Interspec-EpiState-NMF/Interspec-EpiState-NMF.sh /Path_to_Interspec-EpiState-NMF/input_files/config.info.txt 2> test_run.log.txt
```


## Output of Interspec-EpiState-NMF pipeline
- The output target signal matrix after JMnorm should be formatted as N-by-(M+1) matrices, where N represents the number of cCREs, and M represents the number of chromatin features. The first column of each matrix contains the cCRE IDs. The signal values in orignal linear scale for each chromatin feature in the cCREs are saved in the 2~M columns.
- Example output target signal matrix after JMnorm can be found in this [Target.JMnorm_sigmat.txt](https://github.com/guanjue/JMnorm/blob/main/docs/TCD8.JMnorm_sigmat.txt).

## Support
For questions or issues, please either create an issue on the GitHub repository or feel free to reach out via the following email addresses: guanjuexiang@gmail.com

## LICENSE
This project is licensed under the GNU GENERAL PUBLIC License (Version >=2.0). See the [LICENSE](https://github.com/guanjue/JMnorm/blob/main/LICENSE) file for details.







