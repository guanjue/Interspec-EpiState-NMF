cd /Path_to_Interspec-EpiState-NMF/input_files/

##########################################
### for Human gene VS Mouse gene analysis
# download the epigenetic state file 'S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state' from this link: https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/
# download the epigenetic state file 'S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state' from this link: https://usevision.org/data/mm10/ideasJointMay2021/
#
# convert space serperated IDEAS epigenetic state file to bed file (select the columns of shared cell types in both human and mouse) 
cat S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4, $5, $6,$7, $12, $16,$17, $18, $20, $30, $32,$33, $34, $42, $44, $46  }' > S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed
cat S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4, $5, $6,$7, $11, $14,$15, $20, $22, $25, $26,$27, $28, $30, $31, $32  }' > S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed
##########################################


##########################################
### for Human genes analysis
# download the epigenetic state file 'S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state' from this link: https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/
#
# convert space serperated IDEAS epigenetic state file to bed file (remove the first column and last column and 40th and 41st column (remove NEU))
tail -n+2 S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state | cut -d ' ' -f 2- | rev | cut -d ' ' -f 2- | rev \
| sed 's/ /\t/g' | cut -f 1-38,41- > S3V2_IDEAS_hg38_r3_withHg38Mm10prior.all_human_cell_types.state.bed
##########################################


##########################################
### for Mouse genes analysis
# download the epigenetic state file 'S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state' from this link: https://usevision.org/data/mm10/ideasJointMay2021/
#
# convert space serperated IDEAS epigenetic state file to bed file (remove the first column and last column and 40th and 41st column (remove NEU))
tail -n+2 S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state | cut -d ' ' -f 2- | rev | cut -d ' ' -f 2- | rev \
| sed 's/ /\t/g' > S3V2_IDEAS_mm10_r3_withHg38Mm10prior.all_mouse_cell_types.state.bed
##########################################

