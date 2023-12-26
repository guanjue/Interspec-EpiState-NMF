- The `GATA1.Gata1.cor.heatmap.png` and `GATA1.Gata1.cor.heatmap.png.cor.mat.txt` are the correlation matrix of the cross-cell-type epigenetic state patterns between the Human target gene and Mouse target gene. Each column represents a genomic bin in Human gene locus. Each row represents a genomic bin in Mouse gene locus.
```
>>> ls /output_folder/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/GATA1.Gata1.cor.heatmap.png*
GATA1.Gata1.cor.heatmap.png		GATA1.Gata1.cor.heatmap.png.cor.mat.txt
```

- The NMF results are saved in the `NMF_reconstruction_GATA1_Gata1` folder
```
>>> ls /output_folder/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/NMF_reconstruction_GATA1_Gata1
GATA1.Gata1.NMFs.num.BIC.pdf				GATA1.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.txt
GATA1.Gata1.cross_factor_cor_mat.pdf			Gata1.NMFs.num.BIC.pdf
GATA1.Gata1.cross_factor_cor_mat.txt			Gata1.cor.heatmap.png.cor.mat.txt.NMFs.mm10.heatmap.png
GATA1.cor.heatmap.png.cor.mat.txt.NMFs.hg38.heatmap.png	Gata1.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.heatmap.png
GATA1.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.heatmap.png	Gata1.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.txt
```

- The heatmaps of NMF decomposition matrices in the two species are saved in `NMF_reconstruction_GATA1_Gata1/GATA1.cor.heatmap.png.cor.mat.txt.NMFs.hg38.heatmap.png` and `NMF_reconstruction_GATA1_Gata1/Gata1.cor.heatmap.png.cor.mat.txt.NMFs.mm10.heatmap.png`. 
```
ls /output_folder/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/NMF_reconstruction_GATA1_Gata1/*.cor.heatmap.png.cor.mat.txt.NMFs.*png 
NMF_reconstruction_GATA1_Gata1/GATA1.cor.heatmap.png.cor.mat.txt.NMFs.hg38.heatmap.png
NMF_reconstruction_GATA1_Gata1/Gata1.cor.heatmap.png.cor.mat.txt.NMFs.mm10.heatmap.png
```

- The binary heatmap of the identified highly correlated regions based NMF decomposition matrices (FDR-based threshold) in the two species are saved in `NMF_reconstruction_GATA1_Gata1/GATA1.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.heatmap.png` and `NMF_reconstruction_GATA1_Gata1/Gata1.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.heatmap.png`. The corresponding matrices `.txt` are saved in `GATA1.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.txt` and `Gata1.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.txt`. Each row represents a genomic bin in Human / Mouse. Each column represents a NMF factor.
```
>>> # Human gene's FDR high correlation calls
>>> head GATA1.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.txt
H_chrX_48760000_48760200	0	0	0	0	0	0
H_chrX_48760200_48760400	0	0	0	0	0	0
H_chrX_48760400_48760600	0	0	0	0	0	0
H_chrX_48760600_48760800	0	0	0	0	0	0
H_chrX_48760800_48761000	0	0	0	0	0	0
H_chrX_48761000_48761200	0	0	0	0	0	0
H_chrX_48761200_48761400	0	0	0	0	0	0
H_chrX_48761400_48761600	0	0	0	0	0	0
H_chrX_48761600_48761800	0	0	0	0	0	0
H_chrX_48761800_48762000	0	0	0	0	0	0
```
```
>>> # Mouse gene's FDR high correlation calls
>>> head Gata1.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.txt 
M_chrX_7915000_7915200	0	0	0	0	0	0
M_chrX_7915200_7915400	0	0	0	0	0	0
M_chrX_7915400_7915600	0	0	0	0	0	0
M_chrX_7915600_7915800	0	0	0	0	0	0
M_chrX_7915800_7916000	0	0	0	0	0	0
M_chrX_7916000_7916200	0	0	0	0	0	0
M_chrX_7916200_7916400	0	0	0	0	0	0
M_chrX_7916400_7916600	0	0	0	0	0	0
M_chrX_7916600_7916800	0	0	0	0	0	0
M_chrX_7916800_7917000	0	0	0	0	0	0
```

- The `hg38.gene.GATA1.matched_ct.state.bed` and `mm10.gene.Gata1.matched_ct.state.bed` contains the cross-cell-type epigenetic state labels at the Human / Mouse gene locus. 
```
>>> # BED file containing data on Human epigenetic states across various cell types at the Human gene locus
>>> head /output_folder/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/hg38.gene.GATA1.matched_ct.state.bed
chrX	48760000	48760200	0	0	0	0	0	0	0	0	0	0	0
chrX	48760200	48760400	0	0	0	0	2	2	0	0	0	0	0
chrX	48760400	48760600	0	0	0	0	2	2	0	0	0	0	0
chrX	48760600	48760800	0	0	0	0	0	2	0	0	0	0	0
chrX	48760800	48761000	0	3	0	0	0	0	0	0	0	0	0
chrX	48761000	48761200	0	0	0	0	0	0	0	0	0	0	0
chrX	48761200	48761400	0	0	0	0	2	2	0	0	0	0	0
chrX	48761400	48761600	0	0	0	0	0	2	0	0	0	0	0
chrX	48761600	48761800	0	0	0	0	2	2	0	0	0	5	0
chrX	48761800	48762000	9	3	0	0	9	2	0	0	0	9	0
```
```
>>> # BED file containing data on Mouse epigenetic states across various cell types at the Mouse gene locus
>>> head /output_folder/test_cormat_NMF_FDR_pipeline_GATA1_Gata1/mm10.gene.Gata1.matched_ct.state.bed 
chrX	7915000	7915200	0	0	0	0	0	0	0	0	0	0	0	0	0
chrX	7915200	7915400	0	0	0	0	0	0	0	0	0	0	0	0	0
chrX	7915400	7915600	0	0	0	0	0	9	0	0	4	0	0	0	0
chrX	7915600	7915800	9	0	9	9	4	9	9	9	9	9	9	9	9
chrX	7915800	7916000	13	9	9	13	13	13	13	13	13	13	9	9	13	13	9
chrX	7916000	7916200	0	0	0	0	0	9	0	9	9	0	0	3	0
chrX	7916200	7916400	0	0	0	0	0	0	0	0	0	0	0	0	0
chrX	7916400	7916600	0	0	0	0	0	0	0	0	0	0	0	0	0
chrX	7916600	7916800	0	0	0	0	0	0	0	0	0	0	0	0	0
chrX	7916800	7917000	0	0	0	0	3	0	0	9	9	0	0	0	0
```




