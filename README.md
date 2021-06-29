# Systematic Evaluation of DNA Sequence Variations on in vivo Transcription Factor Binding Affinity

## Prerequisite:
The following R packages are required for executing the main code file:   
`seqinr`
`gkmSVM`
`Biostrings`
`BSgenome.Hsapiens.UCSC.hg19`

All the other dependent packages are listed at the beginning of each code file.

## Original Datasets:

(1) "bed_file": all peak files used as positive sets for training gkmSVM. The details can be found in Supplementary Table S2.

(2) "PWM": Source position frequency matrix for target TFs. Most files are from JASPAR, Factorbook and Xu, et al.(2015). All PFMs are transfered to position probability matrix ahead of analysis, thus there are 18 files ended with ".ppm" for use in this folder. The details can be found in supplementary materials.

(3) "SNPs_db153Common": All SNPs found on 23 human chromosomes (from UCSC Genome Browser)

(4) "Disease_SNP": disease related SNPs, selected from all SNPs (based on PheGeni)


## Code Details:
Functions are specified in each part:   

Part 1: Dateset preparation     

Part 3: Prepare datasets for TFs/ Generate control sequences/ Proceed score matching   

Part 4: Threshold determination

Part 5: Position-specific Summary

Part 7: Measuring impact of SNV on TF binding strength (WT and alternative motif preparation)  

Part 9: Exploring potential association between TFs and complex diseases

Part 10: Supplementary materials (Figures & Tables)

Part 11: Additional performance comparison study for gkmSVM and traditional PWM approach.
