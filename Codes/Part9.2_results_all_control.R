library(dplyr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library("writexl")


# read all the data from PheGeni for every disease

wd <- "/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/"
Disease <- c("AD", "Asthma", "Breast neoplasms", "Cardiovascular diseases", "Child Development Disorders Pervasive", "Colorectal Neoplasms", "Crohn diseases", "Lung neoplasms", "Obesity", "Psoriasis", "Type 2 diabetes")



for(i in 1:length(Disease)){
  Di <- Disease[i]
  SNPs <- read.delim(paste0(wd,Di,"/SNPs.csv"), header = T)
  SNPs <- SNPs[!duplicated(SNPs$SNP.rs),]

  save(SNPs,file = paste0(wd,Di,"/SNPs.rda"))

  SNP_lower_region <- SNPs$Location - 5000 + 1000000
  SNP_upper_region <- SNPs$Location + 5000 + 1000000
  SNP_region <- data.frame(SNP = SNPs$SNP.rs, Chr = SNPs$Chromosome, lower = SNP_lower_region, upper = SNP_upper_region)

  save(SNP_region,file = paste0(wd, Di, "/control/SNP_region.rda"))
  
  SNP_region2 <- data.frame(Chr = paste0("chr",SNPs$Chromosome), lower = SNP_lower_region, upper = SNP_upper_region)
  if(i == 1){
    SNP_region2 <- SNP_region2[-58,]
    SNP_region2 <- SNP_region2[-246,]
  } else if(i == 2){
    SNP_region2 <- SNP_region2[-106,]
    SNP_region2 <- SNP_region2[-290,]
    SNP_region2 <- SNP_region2[-290,]
    SNP_region2 <- SNP_region2[-290,]
  } else if(i == 4){
    SNP_region2 <- SNP_region2[-44,]
  } else if(i == 5){
    SNP_region2 <- SNP_region2[-136,]
  } else if(i == 6){
    SNP_region2 <- SNP_region2[-24,]
    SNP_region2 <- SNP_region2[-27,]
    SNP_region2 <- SNP_region2[-33,]
  } else if(i == 7){
    SNP_region2 <- SNP_region2[-88,]
    SNP_region2 <- SNP_region2[-109,]
    SNP_region2 <- SNP_region2[-131,]
    SNP_region2 <- SNP_region2[-291,]
    SNP_region2 <- SNP_region2[-395,]
  } else if(i == 9){
    SNP_region2 <- SNP_region2[-74,]
  } else if(i == 10){
    SNP_region2 <- SNP_region2[-284,]
  } 
  write.table(SNP_region2, paste0(wd, Di, "/control/SNP_region2.txt"), append = F, row.names = F, col.names = F)

}

count_mutate <- function(x){length(unlist(strsplit(as.character(x),"")))}

ref.mutate = function(alle1){
  temp1 = (unlist(strsplit(as.character(alle1),",")))
  return(temp1)
}


for(i in 1:length(Disease)){
  Di <- Disease[i]
  SNP_10kb <- read.delim(paste0(wd,Di,"/control/SNP-10kb.txt"),header = T)
  save(SNP_10kb, file = paste0(wd,Di,"/control/SNP_10kb.rda"))
  
  
  mutate.count = mapply(count_mutate, x = SNP_10kb$ref)
  
  # keep mutate count = 1
  SNP_10kb_mut2 = SNP_10kb[mutate.count == 1,]
  
  # keep those alts with only 1bp
  mutate.count = mapply(count_mutate, x = SNP_10kb_mut2$alts)
  SNP_10kb_mut2 = SNP_10kb_mut2[mutate.count == 2,]
  SNP_10kb_mut2 <- SNP_10kb_mut2[!duplicated(SNP_10kb_mut2$name),]
  
  ref_mut = t(mapply(ref.mutate,alle1 = SNP_10kb_mut2$alts))
  
  SNP_10kb_mut2 = cbind(SNP_10kb_mut2[,1:4],ref1 = as.character(SNP_10kb_mut2$ref),ref2=as.character(ref_mut[1,]))
  save(SNP_10kb_mut2, file = paste0(wd, Di, "/control/SNP_10kb_mut2.rda"))
  
  SNP_10kb_mut2_23 = SNP_10kb_mut2[which(SNP_10kb_mut2$X.chrom %in% paste0("chr",c(1:22,"X"))),]
  save(SNP_10kb_mut2_23, file = paste0(wd, Di, "/control/SNP_10kb_mut2_23.rda"))
  
}




# extracting corresponding sequence from dataset
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

#chromosomes of interest
my_chr <- c(1:22,'X')
my_chr <- gsub(pattern="^", replacement='chr', my_chr)

for(j in 1:length(Disease)){
  Di <- Disease[j]
  SNP_10kb_mut2_23 <- get(load(paste0(wd, Di, "/control/SNP_10kb_mut2_23.rda")))
  
  lab.GRange = GRanges(seqnames = Rle(SNP_10kb_mut2_23$X.chrom), ranges = IRanges(start = SNP_10kb_mut2_23$chromEnd-9, end = SNP_10kb_mut2_23$chromEnd+9))
  
  ref = Hsapiens
  target.snp = list()
  
  for (i in 1:23){
    Chr = paste0("chr",c(1:22,"X"))
    
    
    if(!my_chr[i] %in% lab.GRange@seqnames){
      next
    }
    
    Seq = Hsapiens[[my_chr[i]]]
    thislab = lab.GRange[lab.GRange@seqnames == Chr[i]]
    # thislab
    # length(Seq)        length(Seq) should be greater than thislab
    Seq.set=DNAStringSet(Seq, start=start(thislab[1:length(thislab)]), end=end(thislab[1:length(thislab)]))
    target.snp[[i]] = as.character(Seq.set)
    target.snp_2 <- data.frame(target_snp = unlist(target.snp[[i]]),
                               chr = SNP_10kb_mut2_23[SNP_10kb_mut2_23$X.chrom == my_chr[i],]$X.chrom,
                               position = SNP_10kb_mut2_23[SNP_10kb_mut2_23$X.chrom == my_chr[i],]$chromEnd, 
                               name = SNP_10kb_mut2_23[SNP_10kb_mut2_23$X.chrom == my_chr[i],]$name)
    write.table(target.snp_2, paste0(wd, Di, "/control/target.snp_2.txt"), append = T, row.names = F,col.names = F)
  }
  
  save(target.snp, file = paste0(wd, Di, "/control/target.snp.rda"))
  
  target.snp_2 <- read.table(paste0(wd,Di,"/control/target.snp_2.txt"), header = F)
  colnames(target.snp_2) <- c("seq", "chr", "position","name")
  save(target.snp_2, file = paste0(wd, Di, "/control/target.snp_2.rda"))
  
}



# input: target.snp, ADSNP_10kb_mut2_23
seqmut <- function(target.snp_2, SNP_10kb_mut2_23){
  seq.ref <- list()
  seq.mut <- list()
  for(i in 1:23){
    if(!Chr[i] %in% target.snp_2$chr){
      next
    }
    
    seq.ref[[i]] = target.snp_2[target.snp_2$chr == Chr[i],]$seq
    alleinfo = SNP_10kb_mut2_23[which(SNP_10kb_mut2_23$X.chrom==Chr[i]),]
    seq.mut[[i]] = rep("N",nrow(alleinfo))
    seq.mut[[i]][which(substr(seq.ref[[i]],10,10)==alleinfo$ref1)] = as.character(alleinfo$ref2[which(substr(seq.ref[[i]],10,10)==alleinfo$ref1)])
    seq.mut[[i]][which(substr(seq.ref[[i]],10,10)==alleinfo$ref2)] = as.character(alleinfo$ref1[which(substr(seq.ref[[i]],10,10)==alleinfo$ref2)])
    print(length(which(seq.mut[[i]]=="N")))
  }
  return(seq.mut)
}


for(i in 1:length(Disease)){
  Di <- Disease[i]
  target.snp <- get(load(paste0(wd, Di, "/control/target.snp.rda")))
  target.snp_2 <- get(load(paste0(wd, Di, "/control/target.snp_2.rda")))
  SNP_10kb_mut2_23 <- get(load(paste0(wd, Di, "/control/SNP_10kb_mut2_23.rda")))
    
  seq.mut <- seqmut(target.snp_2, SNP_10kb_mut2_23)
  seq.ref = unlist(target.snp)
  seq.snp = paste0(substr(seq.ref,1,9),unlist(seq.mut),substr(seq.ref,11,19))
  
  seq.ref = data.frame(seq = seq.ref, chr_position = paste(target.snp_2$chr,target.snp_2$position), name = target.snp_2$name)
  seq.snp = data.frame(seq = seq.snp, chr_position = paste(target.snp_2$chr,target.snp_2$position), name = target.snp_2$name)
  
  save(seq.ref,file = paste0(wd, Di, "/control/seq.ref.rda"))
  save(seq.snp,file = paste0(wd, Di, "/control/seq.snp.rda"))
}






# gkmSVM

de = as.matrix(d)

seq10_deltasum = function(x){
  xlist = mapply(substr, x, 1:10, 10:19, USE.NAMES=FALSE)
  xlistr= chartr("ATGC","TACG",xlist)
  
  m1 = d[match(xlist,de[,1]),2]
  m2 = d[match(xlistr,de[,1]),2]
  m1[which(is.na(m1))]=m2[which(is.na(m1))]
  sum = sum(m1)
  return(sum)
}



# threshold of random seq

get.quantile <- function(ii){
  TF = name[ii]
  load(paste0(wd,TF,"/rand.gm.rda"))
  
  gm12878.left  = quantile(c(as.vector(rand.gm)),probs=0.025,na.rm = T)
  gm12878.right = quantile(c(as.vector(rand.gm)),probs=0.975,na.rm = T)
  
  thre = data.frame(TF,gm12878.left,gm12878.right)
  write.table(thre, paste0(wd,"quantile.new.txt"), append = T, row.names = F,col.names = F)
}

for(ii in 1:length(name)){
  get.quantile(ii)
}
quan = read.table(paste0(wd,"quantile.new.txt"), header = F)
colnames(quan) = c("TF", "left", "right")
write.csv(quan, paste0(wd,"quantile.new.csv"), row.names = F)


# for all motifs, get the results
name = c("BCL11A", "CTCF", "EGR1", "GABPA", "JUND", "JUN", "MAX","NANOG", "POU5F1", "RAD21", "REST","RFX5", "SIX5", "SRF","STAT1", "TCF12", "USF1","USF2","YY1")


# get deltaSVM for every overlapped SNP and find significant SNPs


results <- function(Disease){
  
  
  for(i in 1:length(Disease)){
    Di = Disease[i]
    SNP_10kb_mut2_23 <- get(load(paste0(wd, Di, "/control/SNP_10kb_mut2_23.rda")))
    seq.ref <- get(load(paste0(wd, Di, "/control/seq.ref.rda")))
    seq.snp <- get(load(paste0(wd, Di, "/control/seq.snp.rda")))
    
    SNP_10kb_mut2_23_new <- SNP_10kb_mut2_23[,-5:-6]
    colnames(SNP_10kb_mut2_23_new) <- c("seqnames", "start", "end", "snps")
    SNP_10kb_mut2_23_new$start <- SNP_10kb_mut2_23_new$end
    SNP_10kb_mut2_23_new <- makeGRangesFromDataFrame(SNP_10kb_mut2_23_new, keep.extra.columns=T)
    
    
    for(j in 1:length(name)){
      TF = name[j]
      motif.TF <- get(load(paste0("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/motif/", "motif.", TF, ".rdata")))
      
      # matching overlapped SNPs
      x3 <- SNP_10kb_mut2_23_new %over% motif.TF
      if(length(which(x3)) == 0){
        next
      }
      
      overlapSNP <- SNP_10kb_mut2_23_new[which(x3),]
      overlapSNP <- as.data.frame(overlapSNP)
      snpname <- SNP_10kb_mut2_23$name[x3]
      
      # find reference and mutate sequence
      
      seqref <- as.character(seq.ref[match(snpname, seq.ref$name),]$seq)
      seqsnp <- as.character(seq.snp[match(snpname, seq.ref$name),]$seq)
      
      
      
      result <- data.frame(snp = snpname, chr = overlapSNP$seqnames, position = overlapSNP$end, Seqref = seqref, Seqmut = seqsnp)
      
      rownames(result) <- NULL
      result <- result[!grepl("N", result$Seqmut),]
      save(result, file = paste0(wd, Di, "/control/result", TF, ".rda"))
      
      
    }
    
  }
  
}



Disease <- c("AD", "Asthma", "Breast neoplasms", "Cardiovascular diseases", "Child Development Disorders Pervasive", "Colorectal Neoplasms", "Crohn diseases", "Lung neoplasms", "Obesity", "Psoriasis", "Type 2 diabetes")
results(Disease)


for(i in 1:length(Disease)){
  Di <- Disease[i]
  for(j in 1:length(name)){
    TF <- name[j]
    f1 <- paste0(wd, Di, "/control/result", TF, ".rda")
    if(!file.exists(f1)){
      next
    }
    
    result <- get(load(paste0(wd, Di, "/control/result", TF, ".rda")))
    d <- get(load(paste0(wd, "motif/", TF, "/d.rda")))
    de <- as.matrix(d)
    
    gkmSVM.ref = mapply(seq10_deltasum,result$Seqref)
    gkmSVM.snp = mapply(seq10_deltasum,result$Seqmut)
    save(gkmSVM.ref, file = paste0(wd,Di, "/control/gkmSVM.ref", TF, ".rda"))
    save(gkmSVM.snp, file = paste0(wd,Di, "/control/gkmSVM.snp", TF, ".rda"))
    
    deltaSVM = gkmSVM.snp - gkmSVM.ref
    
    result <- cbind(result, deltaSVM)
    result$sig_SNPs <- ifelse(result$deltaSVM < quan$left[j] | result$deltaSVM > quan$right[j], 1 , 0)
    save(result, file = paste0(wd, Di, "/control/result", TF, ".rda"))
    write.table(result, paste0(wd, Di, "/control/result_", TF, ".txt"), append = F, row.names = F,col.names = F)
  }
}


# remove duplicates in the result for each motif
for(i in 1:length(Disease)){
  Di <- Disease[i]
  for(j in 1:length(name)){
    TF <- name[j]
    f1 <- paste0(wd, Di, "/control/result", TF, ".rda")
    if(!file.exists(f1)){
      next
    }
    result_motifs <- get(load(paste0(wd, Di, "/control/result", TF, ".rda")))
    result_motifs <- result_motifs[order(result_motifs$deltaSVM, decreasing = T),]
    result_motifs <- result_motifs[!duplicated(result_motifs$snp),]
    save(result_motifs, file = paste0(wd, Di, "/control/result_motifs", TF, ".rda"))
  }
}



for(i in 1:length(Disease)){
  Di <- Disease[i]
  final_result = read.table(paste0(wd, Di, "/control/result.txt"), header = F)
  colnames(final_result) = c("snp", "chr", "position", "Seqref", "Seqmut", "DeltaSVM", "significant SNP")
  write.csv(final_result, paste0(wd, Di, "/control/result.csv"), row.names = F)
  write_xlsx(final_result,path = paste0(wd, Di, "/control/result.xlsx"))
  
}



for (i in 1:length(Disease)) {
  Di <- Disease[i]
  
  for(j in 1:length(name)){
    TF <- name[j]
    wd2 <- paste0(wd, Di, "/control/result", TF, ".rda")
    if(!file.exists(wd2)) {
      next
    }
    
    result.TF <- get(load(paste0(wd, Di, "/control/result", TF, ".rda")))
    
    result_1st <- cbind(nrow(result.TF),result.TF[1,])
    colnames(result_1st) <- c("nrow", colnames(result.TF))
    write.table(result_1st, paste0(wd, Di, "/control/result_1st.txt"), append = T, row.names = F,col.names = F)
  }
  
}





resultBCL11A <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/BCL11A/resultBCL11A.rda"))
resultCTCF <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/CTCF/resultCTCF.rda"))
resultEGR1 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/EGR1/resultEGR1.rda"))
resultGABPA <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/GABPA/resultGABPA.rda"))
resultJUND <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/JUND/resultJUND.rda"))
resultJUN <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/JUN/resultJUN.rda"))
resultMAX <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/MAX/resultMAX.rda"))
resultNANOG <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/NANOG/resultNANOG.rda"))
resultPOU5F1 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/POU5F1/resultPOU5F1.rda"))
resultRAD21 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/RAD21/resultRAD21.rda"))
resultRFX5 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/RFX5/resultRFX5.rda"))
resultSIX5 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/SIX5/resultSIX5.rda"))
resultSRF <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/SRF/resultSRF.rda"))
resultSTAT1 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/STAT1/resultSTAT1.rda"))
resultTCF12 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/TCF12/resultTCF12.rda"))
resultUSF1 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/USF1/resultUSF1.rda"))
resultUSF2 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/USF2/resultUSF2.rda"))
resultYY1 <- get(load("/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/AD SNPs/YY1/resultYY1.rda"))

nrow(resultBCL11A) #543
nrow(resultCTCF) #14
nrow(resultEGR1) #218
nrow(resultGABPA) #208
nrow(resultJUND) #418
nrow(resultJUN) #229
nrow(resultMAX) #583
nrow(resultNANOG) #543
nrow(resultPOU5F1) #1384
nrow(resultRAD21) #94
nrow(resultRFX5) #171
nrow(resultSIX5) #45
nrow(resultSRF) #88
nrow(resultSTAT1) #275
nrow(resultTCF12) #733
nrow(resultUSF1) #235
nrow(resultUSF2)  #222
nrow(resultYY1) #200

resultBCL11A[1,]   # rs112727566 chr19 44893991 TGTCAGCTTGACATTGTTG TGTCAGCTTAACATTGTTG 0.186655        0
resultCTCF[1,]     # rs188388424 chr11 60259290 GCACCACCATACCCAGCTA GCACCACCACACCCAGCTA -0.970647        0
resultEGR1[1,]     # rs34662182 chr11 86159223 CCTGTGAGCTGCTTCCTCC CCTGTGAGCAGCTTCCTCC -0.546785        0
resultGABPA[1,]    # rs59953894 chr19 44890841 AAAAAAAAAAAAAAACACT AAAAAAAAACAAAAACACT -2.315062        0
resultJUND[1,]     # rs148795043 chr19 44890674 TCAAGGTTGTGGTGAGCCA TCAAGGTTGCGGTGAGCCA -0.164122        0
resultJUN[1,]      # rs2722664 chr19 44909990 GGTTGCAGTGAGTTGAAAC GGTTGCAGTAAGTTGAAAC 0.998251        0
resultMAX[1,]      # rs2722640 chr19 44892622 AAAATGAGACGCCTGTTTC AAAATGAGATGCCTGTTTC 0.599994        0
resultNANOG[1,]    # rs112727566 chr19 44893991 TGTCAGCTTGACATTGTTG TGTCAGCTTAACATTGTTG 0.186655        0
resultPOU5F1[1,]   # rs12978706 chr19 44893523 GCAAGATGCAGTCTGGGAG GCAAGATGCTGTCTGGGAG -0.933397        0
resultRAD21[1,]    # rs116761273 chr1 154451723 ATATTTTCTAATTTTTCCA ATATTTTCTGATTTTTCCA  0.35978        0
resultRFX5[1,]     # rs2490256 chr1 207521901 CCTGCCAGCTCTGCGTGAA CCTGCCAGCGCTGCGTGAA 1.491266        0
resultSIX5[1,]     # rs2722755 chr19 44882699 GCTGGGTGATGGTTTTATA GCTGGGTGACGGTTTTATA -0.936297        0
resultSRF[1,]      # rs4147922 chr19  1062598 ATTATGTCACTGTTCACTA ATTATGTCAGTGTTCACTA -1.230701        0
resultSTAT1[1,]    # rs117829793 chr19 44923735 AGATAAAAACTATCAGTGT AGATAAAAAATATCAGTGT 0.944798        0
resultTCF12[1,]    # rs8109490 chr19 44893577 ACGGAATCTCGCTCTGTCA ACGGAATCTTGCTCTGTCA -1.195169        0
resultUSF1[1,]     # rs139135897 chr19 44895375 GTGCTGAAAATTTTTTGTC GTGCTGAAATTTTTTTGTC -1.416925        0
resultUSF2[1,]     # rs112275114 chr19 44895424 CCAAAAACCGTATATTCCT CCAAAAACCATATATTCCT -1.123084        0
resultYY1[1,]      # rs112642778 chr2 127133869 CACATAGCACAAGGTCTGG CACATAGCATAAGGTCTGG -0.441219        0





