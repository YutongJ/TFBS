options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)
arguments = matrix(unlist(strsplit(args,"=")),ncol=2,byrow = T)

for (args_i in 1:length(args)) {
  assign(arguments[args_i,1],as.numeric(arguments[args_i,2]))
}




name = c("BCL11A", "CTCF", "EGR1", "GABPA", "JUND", "JUN", "MAX","NANOG", "POU5F1", "RAD21", "REST","RFX5", "SIX5", "SRF","STAT1", "TCF12", "USF1","USF2","YY1")
wd <- "/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/"
Disease <- c("AD", "Asthma", "Breast neoplasms", "Cardiovascular diseases", "Child Development Disorders Pervasive", "Colorectal Neoplasms", "Crohn diseases", "Lung neoplasms", "Obesity", "Psoriasis", "Type 2 diabetes")

results_not_overlap <- function(Disease){
  
  
  for(i in 1:length(Disease)){
    Di = Disease[i]
    SNP_10kb_mut2_23 <- get(load(paste0(wd, Di, "/SNP_10kb_mut2_23.rda")))
    seq.ref <- get(load(paste0(wd, Di, "/seq.ref.rda")))
    seq.snp <- get(load(paste0(wd, Di, "/seq.snp.rda")))
    
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
        not_overlapSNP <- SNP_10kb_mut2_23_new
      } else{
        not_overlapSNP <- SNP_10kb_mut2_23_new[-which(x3),]
      }
      
      not_overlapSNP <- as.data.frame(not_overlapSNP)
      snpname <- not_overlapSNP$snps
      
      # find reference and mutate sequence
      
      seqref <- as.character(seq.ref[match(snpname, seq.ref$name),]$seq)
      seqsnp <- as.character(seq.snp[match(snpname, seq.ref$name),]$seq)
      
      
      
      result_not_overlap <- data.frame(snp = snpname, chr = not_overlapSNP$seqnames, position = not_overlapSNP$end, Seqref = seqref, Seqmut = seqsnp)
      
      rownames(result_not_overlap) <- NULL
      result_not_overlap <- result_not_overlap[!grepl("N", result_not_overlap$Seqmut),]
      save(result_not_overlap, file = paste0(wd, Di, "/result_not_overlap_", TF, ".rda"))
      write.csv(result_not_overlap, file = paste0(wd, Di, "/result_not_overlap_", TF, ".csv"))
      
    }
    
  }
  
}

Disease <- c("AD", "Asthma", "Breast neoplasms", "Cardiovascular diseases", "Child Development Disorders Pervasive", "Colorectal Neoplasms", "Crohn diseases", "Lung neoplasms", "Obesity", "Psoriasis", "Type 2 diabetes")
results_not_overlap(Disease)


# calculate deltaSVM
# 
# for(i in 1:length(Disease)){
#   Di <- Disease[i]
#   for(j in 1:length(name)){
#     TF <- name[j]
#    
#     
#     result_not_overlap <- get(load(paste0(wd, Di, "/result_not_overlap_", TF, ".rda")))
#     d <- get(load(paste0(wd, "motif/", TF, "/d.rda")))
#     de <- as.matrix(d)
#     
#     gkmSVM.ref.not.overlap = mapply(seq10_deltasum,result_not_overlap$Seqref)
#     gkmSVM.snp.not.overlap = mapply(seq10_deltasum,result_not_overlap$Seqmut)
#     save(gkmSVM.ref.not.overlap, file = paste0(wd,Di, "/gkmSVM.ref.not.overlap_", TF, ".rda"))
#     save(gkmSVM.snp.not.overlap, file = paste0(wd,Di, "/gkmSVM.snp.not.overlap_", TF, ".rda"))
#     
#     deltaSVM = gkmSVM.snp.not.overlap - gkmSVM.ref.not.overlap
# 
#     result_not_overlap <- cbind(result_not_overlap[,1:5], deltaSVM)
#     result_not_overlap$sig_SNPs <- ifelse(deltaSVM < quan$left[j] | deltaSVM > quan$right[j], 1 , 0)
#     save(result_not_overlap, file = paste0(wd, Di, "/result_not_overlap_", TF, ".rda"))
#     write.table(result_not_overlap, paste0(wd, Di, "/result_not_overlap_", TF, ".txt"), append = F, row.names = F,col.names = F)
#   }
# }



for(i in 1:length(Disease)){
  Di <- Disease[i]
  for(j in 1:length(name)){
    TF <- name[j]
    deltaSVM <- get(load(paste0(wd, Di, "/deltaSVM.not.overlap_", TF, ".rda")))
    result_not_overlap <- get(load(paste0(wd, Di, "/result_not_overlap_", TF, ".rda")))
    result_not_overlap <- cbind(result_not_overlap[,1:5], deltaSVM)
    result_not_overlap$sig_SNPs <- ifelse(deltaSVM < quan$left[1] | deltaSVM > quan$right[1], 1 , 0)
    save(result_not_overlap, file = paste0(wd, Di, "/result_not_overlap_sig", TF, ".rda"))
    
  }
  
}


for(i in 1:length(Disease)){
  Di <- Disease[i]
  SNP_10kb_mut2_23 <- get(load(paste0(wd, Di, "/SNP_10kb_mut2_23.rda")))
  
  for(j in 1:length(name)){
    TF <- name[j]
    f1 <- paste0(wd, Di, "/result", TF, ".rda")
    if(!file.exists(f1)){
      next
    }
    result <- get(load(paste0(wd, Di, "/result", TF, ".rda")))
    result_not_overlap_sig <- get(load(paste0(wd, Di, "/result_not_overlap_sig", TF, ".rda")))
    
    FP_FN <- data.frame(total_SNP = nrow(SNP_10kb_mut2_23),
                        overlap_with_PWM_motif_sites = nrow(result),
                        overlapped_and_significant_deltaSVM = sum(result$sig_SNPs),
                        False_postive = (nrow(result)-sum(result$sig_SNPs))/nrow(result),
                        not_overlap_with_PWM_motif_sites = nrow(SNP_10kb_mut2_23) - nrow(result),
                        not_overlapped_but_sigificant_deltaSVM = sum(result_not_overlap_sig$sig_SNPs),
                        False_negative = sum(result_not_overlap_sig$sig_SNPs)/(nrow(SNP_10kb_mut2_23) - nrow(result)))
    write.table(FP_FN, file = paste0(wd, Di, "/FP_FN.txt"), append = T, row.names = F,col.names = F)
    
  }
  
  FP_FN <- read.table(paste0(wd, Di, "/FP_FN.txt"), header = F)
  if(nrow(FP_FN) == 19){
    FP_FN_result <- cbind(name, FP_FN)
  } else if(nrow(FP_FN) == 18){
    FP_FN_result <- cbind(name[-11], FP_FN)
  }
  colnames(FP_FN_result) <- c("motif", 
                              "total_SNP", 
                              "overlap_with_PWM_motif_sites", 
                              "overlapped_and_significant_deltaSVM", 
                              "False_postive",
                              "not_overlap_with_PWM_motif_sites",
                              "not_overlapped_but_sigificant_deltaSVM",
                              "False_negative")
  write.csv(FP_FN_result, file = paste0(wd, Di, "/FP_FN.csv"))
  
}


require(openxlsx)
tables = list()
for(i in 1:length(Disease)){
  Di <- Disease[i]
  tmp <- read.csv(paste0(wd, Di, "/FP_FN.csv"), header = T)
  tables[[i]] <- tmp
}

names(tables) <- Disease
write.xlsx(tables, file = paste0(wd, "FP_FN.xlsx"), colWidths = rep("auto", length(Disease)))


tables_FP <- data.frame(motif = tables[[1]][-11,]$motif, AD = tables[[1]][-11,]$False_postive, Asthma = tables[[2]][-11,]$False_postive, Breast_neoplasms = tables[[3]][-11,]$False_postive, Cardiovascular = tables[[4]]$False_postive,
                       CDDP = tables[[5]]$False_postive, Colorectal = tables[[6]][-11,]$False_postive, Crohn = tables[[7]][-11,]$False_postive, Lung_neoplasms = tables[[8]]$False_postive, Obesity = tables[[9]]$False_postive,
                       Psoriasis = tables[[10]][-11,]$False_postive, Type_2_diabetes = tables[[11]]$False_postive)

tables_FN <- data.frame(motif = tables[[1]][-11,]$motif, AD = tables[[1]][-11,]$False_negative, Asthma = tables[[2]][-11,]$False_negative, Breast_neoplasms = tables[[3]][-11,]$False_negative, Cardiovascular = tables[[4]]$False_negative,
                            CDDP = tables[[5]]$False_negative, Colorectal = tables[[6]][-11,]$False_negative, Crohn = tables[[7]][-11,]$False_negative, Lung_neoplasms = tables[[8]]$False_negative, Obesity = tables[[9]]$False_negative,
                            Psoriasis = tables[[10]][-11,]$False_negative, Type2_diabetes = tables[[11]]$False_negative)

tables_FP <- melt(tables_FP)
tables_FN <- melt(tables_FN)

sum_FP <- rep(NA,19)
for(i in 1:length(name)){
  sum_FP[i] <- sum(tables_FP[which(tables_FP$motif == name[i]),]$value)
}

sum_FP <- sum_FP[-11]
a <- data.frame(name = name[-11],sum_FP)
rank <- a[order(a$sum_FP),]

rank2 <- data.frame()
for(i in 1:18){
  rank2 <- rbind(rank2, tables_FP[which(tables_FP$motif == unfactor(rank$name[i])),])
}
rank2 <- rank2[nrow(rank2):1,]

b <- c()
for(i in 1:18){
  b <- c(b, rep(i,11))
}
b <- factor(b,levels = c(1:18), labels = unique(rank2$motif))
rank2$b <- b



library(ggplot2)
library(reshape)
library(lattice)
library(gridExtra)
library(grid)
heatmap_FP <- ggplot(rank2, aes(x = b, y = variable)) + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "steelblue")
heatmap_FP <- heatmap_FP + labs(title = "a")
heatmap_FN <- ggplot(tables_FN, aes(x = motif, y = variable)) + geom_tile(aes(fill = value))
heatmap_FN <- heatmap_FN + labs(title = "b")

panels <- list()
panels[[1]] <- heatmap_FP 
panels[[2]] <- heatmap_FN

pdf(file = "/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/heatmap_FP_FN.pdf", width = 12)
grid.arrange(panels[[1]], panels[[2]], nrow = 2)
dev.off()

pdf(file = "/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/heatmap_FP2.pdf", width = 12)
grid.arrange(panels[[1]],nrow = 1)
dev.off()







