# Systematic Evaluation of DNA Sequence Variations on in vivo Transcription Factor Binding Affinity

## Prerequisite:
The following R packages are required for executing the main code file:   
`seqinr`
`gkmSVM`
`Biostrings`
`BSgenome.Hsapiens.UCSC.hg19`

All the other dependent packages are listed at the beginning of each code file.

## Original Datasets:

(1) "bed_file": all peak files used as positive sets for training gkmSVM 

(2) "PWM": position weight matrix for target TFs

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


