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



### aply threshold to select not-overlapped but significant SNP ###
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


### calculate false positive and false negative ###
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



##### plot heatmap V2 06/28 #######
for(i in 1:length(Disease)){
  Di <- Disease[i]
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


tables_FP <- data.frame(motif = tables[[1]][-11,]$motif, AD = tables[[1]][-11,]$False_postive, Asthma = tables[[2]][-11,]$False_postive, Breast_neoplasms = tables[[3]][-11,]$False_postive, Cardiovascular = tables[[4]]$False_postive,
                        CDDP = tables[[5]]$False_postive, Colorectal = tables[[6]][-11,]$False_postive, Crohn = tables[[7]][-11,]$False_postive, Lung_neoplasms = tables[[8]]$False_postive, Obesity = tables[[9]]$False_postive,
                        Psoriasis = tables[[10]][-11,]$False_postive, Type_2_diabetes = tables[[11]]$False_postive)

tables_FN <- data.frame(motif = tables[[1]][-11,]$motif, AD = tables[[1]][-11,]$False_negative, Asthma = tables[[2]][-11,]$False_negative, Breast_neoplasms = tables[[3]][-11,]$False_negative, Cardiovascular = tables[[4]]$False_negative,
                        CDDP = tables[[5]]$False_negative, Colorectal = tables[[6]][-11,]$False_negative, Crohn = tables[[7]][-11,]$False_negative, Lung_neoplasms = tables[[8]]$False_negative, Obesity = tables[[9]]$False_negative,
                        Psoriasis = tables[[10]][-11,]$False_negative, Type2_diabetes = tables[[11]]$False_negative)


tables_FP <- melt(tables_FP)
tables_FN <- melt(tables_FN)


### sort according to the value ###
### have the darker column listed first ###

sum_FP <- rep(NA,19)
for(i in 1:length(name)){
  sum_FP[i] <- sum(tables_FP[which(tables_FP$motif == name[i]),]$value)
}

sum_FP <- sum_FP[-11]
a <- data.frame(name = name[-11],sum_FP)
rank <- a[order(a$sum_FP),]

rank2 <- data.frame()
for(i in 1:18){
  rank2 <- rbind(rank2, tables_FP[which(tables_FP$motif == rank$name[i]),])
}
rank2 <- rank2[nrow(rank2):1,]

b <- c()
for(i in 1:18){
  b <- c(b, rep(i,11))
}
b <- factor(b,levels = c(1:18), labels = unique(rank2$motif))
rank2$b <- b
colnames(rank2) <- c("motif", "Diseases", "Value", "TFs")

rank2 <- rank2[rev(order(rank2$Diseases)),]

rank2$Diseases <- as.character(rank2$Diseases)

rank2$Diseases <- c(rep("Type 2 Diabetes", 18),
                    rep("Psoriasis",18),
                    rep("Obesity",18),
                    rep("Lung Cancer", 18),
                    rep("Crohn's Disease", 18),
                    rep("Colorectal Cancer", 18),
                    rep("CDDP", 18),
                    rep("Cardiovascular Diseases", 18),
                    rep("Breast Cancer", 18),
                    rep("Asthma", 18),
                    rep("AD", 18))



rank2$Diseases <- factor(rank2$Diseases , levels = c("Type 2 Diabetes","Psoriasis","Obesity","Lung Cancer","Crohn's Disease","Colorectal Cancer",
                                                     "CDDP","Cardiovascular Diseases","Breast Cancer","Asthma","AD"))


library(ggplot2)
library(reshape)
library(lattice)
library(gridExtra)
library(grid)
heatmap_FP <- ggplot(rank2, aes(x = TFs, y = Diseases)) + geom_tile(aes(fill = Value)) + scale_fill_gradient(low = "white", high = "steelblue")


heatmap_FP <- heatmap_FP  + 
  theme(axis.text.x=element_text(angle = 30, size = 12, vjust = 1, hjust = 1, face = "bold", color = "black")) + 
  theme(axis.text.y=element_text(size = 12, vjust = 0, hjust = 1, face = "bold", color = "black")) +
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "black")) +
  theme(axis.title.y = element_text(size = 14, face = "bold", color = "black")) +
  theme(legend.title = element_text(color = "black", size = 12, face = "bold"), legend.text = element_text(color = "black", size = 12, face = "bold"))


heatmap_FN <- ggplot(tables_FN, aes(x = motif, y = Disease)) + geom_tile(aes(fill = value))
heatmap_FN <- heatmap_FN + labs(title = "b")

panels <- list()
panels[[1]] <- heatmap_FP 
panels[[2]] <- heatmap_FN

pdf(file = "/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/heatmap_FP_FN.pdf", width = 12)
grid.arrange(panels[[1]], panels[[2]], nrow = 2)
dev.off()

pdf(file = "/Users/jiangjiahui/Desktop/emory/Thesis/JJH part/heatmap_FP3.pdf", width = 12)
grid.arrange(panels[[1]],nrow = 1)
dev.off()

