
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

- The details about the input files can be found [Here](https://github.com/guanjue/public_log_descriptions/tree/main/Interspec-EpiState-NMF/Input_file_description):



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







