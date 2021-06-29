library(BSgenome.Hsapiens.UCSC.hg19)
# library(gkmSVM)
library(seqinr)

dir <- "H:/10-mers/20210521_test_CTCF_on_different_ChIPseq/"
setwd(dir)

# read in CTCF from a new dataset ENCFF085HTY.bed
bed <- as.data.frame(read.table(paste0(dir, "ENCFF085HTY.bed"),
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

CTCF_new_bed <- bed[,1:3]
CTCF_new_bed <- CTCF_new_bed[order(CTCF_new_bed[,1]),]
write.table(CTCF_new_bed,
            file=paste0(dir, "ENCFF085HTYtest.bed"), quote=F, sep="\t", row.names=F, col.names=F)

# generate Nullseq
inputBedFN = paste0("ENCFF085HTYtest.bed")
genNullSeqs(inputBedFN)
inputBedFN = CTCF_new_bed
genNullSeqs("H:/10-mers/20210521_test_CTCF_on_different_ChIPseq/ENCFF085HTYtest.bed",
            genomeVersion = 'hg19')

#----------------
# prepare neg/pos set
#----------------
set_prep <- function(file){
  # read neg Set
  negSet <- read.fasta(file)
  neg_list0 <- strsplit(names(negSet),"_")
  # remove fix term
  len <- do.call(c, lapply(neg_list0, length))
  # form bed format
  neg_bed <- do.call(rbind, neg_list0[which(len == 5)])
  # generate list by chr
  neg_list <- as.data.frame(neg_bed[,1:3])
  colnames(neg_list) <- c("seqnames", "start", "end")

  # split the whole set into list by Chr
  list_neg <- split(neg_list, f = as.factor(neg_list$seqnames)) # Split data
  list_neg <- list_neg[which(names(list_neg)%in% paste0("chr",c(1:22,"X")))]
  return(list_neg)
}

# neg Set
list_neg <- set_prep("negSet.fa")
# pos Set
list_pos <- set_prep("posSet.fa")

save(list_neg, file = "list_neg.rda")
save(list_pos, file = "list_pos.rda")

# prepare CTCF pos set
CTCF_new_bed <- as.data.frame(read.table(paste0(dir, "ENCFF085HTYtest.bed"),
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(CTCF_new_bed) <- c("seqnames", "start", "end")
# split the whole set into list by Chr
data_list <- split(CTCF_new_bed, f = as.factor(CTCF_new_bed$seqnames)) # Split data



#----------------
# two functions  - find out overlapped number
#----------------
#--------------------------------------------------
# calculate top-weight-rank tenmer overlapping with each chr
chr_overlap <- function(tenmer, x){
  # make Granges
  psn <- makeGRangesFromDataFrame(x, ignore.strand=T)
  # get seq from granges
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  # turn to characters
  seq_chr <- as.character(seq)
  # find overlapped peaks
  seq_chr_overlap <- grep(tenmer, seq_chr)
  num_chr_overlap <- length(seq_chr_overlap)
  return(num_chr_overlap)
}
# test
# test <- chr_overlap("CCACTAGGGG", data_list[[2]])


#--------------------------------------------------
# input:
# chrlist: list of 23 chrs of peak regions
# tenmer_seq: top-ranked tenmer in certain motif
# output:
# the number of certain tenmer showing up in chrs
#--------------------------------------------------
tenmer_num <- function(chrlist = data_list, tenmer_seq){
  # calculate number for each chr
  num_seq <- mapply(x = chrlist, chr_overlap,MoreArgs = list(tenmer = tenmer_seq))
  # sum up all overlapped
  val <- sum(num_seq)
  
  return(val)
}

# test
# tenmer_num(tenmer_seq = "CCACTAGGGG", chrlist = list_pos)




#----------------
# prepare top-ranked tenmers
#----------------
# load d matrix for CTCF
load("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/CTCF/Part1/d.rda")

# order the sequence based on weights
seq_list <- d[order(d$V2,decreasing = TRUE),]
# we are interested in 1000 top-ranked tenmers
num <- 1000
seq_list_sub <- as.character(seq_list[1:num,1])

tab_tenmer <- cbind(tenmer=seq_list_sub,
                    weights=seq_list[1:num,2],
                    pos_count=NA,
                    neg_count=NA,
                    chisq.test=NA)

# sum(do.call(c, lapply(list_pos, nrow))) # 40790
# sum(do.call(c, lapply(list_neg, nrow))) # 31723

for (i in seq_along(seq_list_sub)){
  pos <- tenmer_num(tenmer_seq = seq_list_sub[i],
                    chrlist = list_pos)
  neg <- tenmer_num(tenmer_seq = seq_list_sub[i],
                    chrlist = list_neg)
  # chi-square test
  mat <- matrix(c(pos, 40790-pos,
                  neg, 31723-neg), ncol=2)
  chitest <- chisq.test(mat)
  tab_tenmer[i,3:5] <- c(pos, neg, chitest$p.value)
}

write.csv(tab_tenmer, file = "CTCF_test_on_new_ChIPseq.csv", row.names = F)

tab_tenmer <- read.csv("CTCF_test_on_new_ChIPseq.csv")

fisher <- function(pos, neg){
  # pos <- x[3];neg <- x[4]
  mat <- matrix(c(pos, 40790-pos,
                  neg, 31723-neg), ncol=2)
  fishertest <- fisher.test(mat)$p.value
  return(fishertest)
}

fisher.test <- mapply(pos=tab_tenmer[,3], neg=tab_tenmer[,4],fisher)

tab_tenmer <- cbind(tab_tenmer, fisher.test)

sum(tab_tenmer$fisher.test<0.05) #819
sum(tab_tenmer$chisq.test<0.05) #809

write.csv(tab_tenmer, file = "CTCF_on_new_ChIPseq.csv", row.names = F)




## ----------------
#  negSet_1000 and posSet_1000
## ----------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)

dir <- "H:/10-mers/20210521_test_CTCF_on_different_ChIPseq/"
setwd(dir)

# original input bed pos file
bed <- as.data.frame(read.table(paste0(dir, "ENCFF085HTY.bed"),
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

# read in negSet (all 33796 reads)
# file = "negSet.fa";number=1000
# file = "posSet.fa";number=1000
set_prep <- function(file, number,lens=350){
  # read neg Set
  negSet <- read.fasta(paste0(dir,"CTCF/",file))
  neg_list0 <- strsplit(names(negSet),"_")
  # remove fix term
  len <- do.call(c, lapply(neg_list0, length))
  # form bed format
  neg_bed <- do.call(rbind, neg_list0[which(len == 5)])
  # generate list by chr
  neg_list <- as.data.frame(neg_bed[,1:3])
  neg_list[,4] <- as.numeric(neg_list[,3]) - as.numeric(neg_list[,2])+1
  colnames(neg_list) <- c("seqnames", "start", "end", "width")
  
  if (file == "negSet.fa"){
    output <- neg_list[1:number,]
  }else{
    bed1 <- bed[which(bed$V3-bed$V2==lens),]
    bed2 <- bed1[nrow(bed1):1,]
    # output <- neg_list[match(bed1[1:1000,2]+1,as.numeric(neg_list$start)),]
    output <- neg_list[match(bed2[1:1000,2]+1,as.numeric(neg_list$start)),]
  }
  
  return(output)
}

# neg Set 1000
list_neg_1000 <- set_prep("negSet.fa",number=1000)

# pos Set 1000
list_pos_1000 <- set_prep("posSet.fa",1000,lens=244)

save(list_neg_1000, file = "list_neg_1000.rda")
save(list_pos_1000, file = "list_pos_1000.rda")



# weights/pwm score
load(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/CTCF/Part1/d.rda"))
de <- as.matrix(d)
dir.pwm <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/"
pwmfn <- paste0(dir.pwm, "CTCF", ".ppm.txt")

# prepare pwm(+) and pwm(-)
pwm = t(as.matrix(read.table(pwmfn)))
pwm = (pwm + 0.001)/colSums(pwm + 0.001)

rownames(pwm) = c("A","C","G","T")
prob1 = pwm
rownames(prob1)=c("A","C","G","T")
pr1=(prob1)
# complimentary one
pr2=pr1[4:1,ncol(pr1):1]
rownames(pr2)=c("A","C","G","T")
rm(pwm,prob1,pwmfn)


# x=list_neg_1000[1:5,]
table_weight_prob = function(x, lens=350){
  # make Granges
  psn <- makeGRangesFromDataFrame(x, ignore.strand=T)
  # get seq from granges
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  # turn to characters
  seq_chr <- as.character(seq)
  
  # t=seq_chr[1]
  match_weight_or_pwm <- function(t){
    seq_vector <- mapply(substr, t, 1:(lens-9), 10:lens, USE.NAMES=FALSE)
    seq_vector_r <- chartr("ATGC","TACG",seq_vector)
    
    # weights
    m1 = d[match(seq_vector,de[,1]),2]
    m2 = d[match(seq_vector_r,d[,1]),2]
    m1[which(is.na(m1))]=m2[which(is.na(m1))]
    val_weight <- m1
    
    # log prob
    m1 = d[match(seq_vector,de[,1]),3]
    m2 = d[match(seq_vector_r,d[,1]),3]
    m1[which(is.na(m1))]=m2[which(is.na(m1))]
    val_prob <- m1
    
    return(list(val_weight, val_prob))
  }
  
  xlist <- mapply(FUN=match_weight_or_pwm, seq_chr, SIMPLIFY = FALSE)
  xlist_weight <- do.call(rbind,lapply(xlist, function(x){x[[1]]}))
  xlist_prob <- do.call(rbind,lapply(xlist, function(x){x[[2]]}))

  return(list(xlist_weight, xlist_prob))
}

# negSet_1000
a = proc.time()
res_neg <- table_weight_prob(x=list_neg_1000[1:1000,])
proc.time() - a

a = proc.time()
res_pos <- table_weight_prob(x=list_pos_1000[1:1000,])
proc.time() - a

# res_rowmean
res_rowmean <- cbind(rep(0:1, each=1000),
                     rbind(do.call(cbind, lapply(res_neg, function(x){rowMeans(x,na.rm = TRUE)})),
                           do.call(cbind, lapply(res_pos, function(x){rowMeans(x,na.rm = TRUE)}))))
colnames(res_rowmean) <- c("group", "gkmSVM", "log_prob")
res_rowmean <- as.data.frame(res_rowmean)

# res_max
res_rowmax <- cbind(rep(0:1, each=1000),
                     rbind(do.call(cbind, lapply(res_neg, function(x){apply(x,1,function(t){max(t, na.rm = TRUE)})})),
                           do.call(cbind, lapply(res_pos, function(x){apply(x,1,function(t){max(t, na.rm = TRUE)})}))))
colnames(res_rowmax) <- c("group", "gkmSVM", "log_prob")
res_rowmax <- as.data.frame(res_rowmax)


# ROC / PRC plots
library(precrec)
library(ggplot2)
library(gridExtra)
library(ROCR)


precrec_obj1 <- evalmod(scores = res_rowmean$log_prob, labels = res_rowmean$group)
precrec_obj2 <- evalmod(scores = res_rowmean$gkmSVM, labels = res_rowmean$group)
precrec::auc(precrec_obj1)
precrec::auc(precrec_obj2)

pdf(file = "ROC.pdf", width = 12, height = 12)
grid.arrange(
autoplot(precrec_obj1, curvetype = "ROC"),
autoplot(precrec_obj1, curvetype = "PRC"),
autoplot(precrec_obj2, curvetype = "ROC"),
autoplot(precrec_obj2, curvetype = "PRC"),nrow=2)
dev.off()

autoplot(precrec_obj1, curvetype = "ROC", show_cb=FALSE)



library(ROCR)


# row means
precrec_obj1 <- evalmod(scores = res_rowmean$log_prob, labels = res_rowmean$group)
precrec_obj2 <- evalmod(scores = res_rowmean$gkmSVM, labels = res_rowmean$group)
precrec::auc(precrec_obj1)
precrec::auc(precrec_obj2)
pred1 <- prediction(res_rowmean$log_prob, res_rowmean$group)
perf11 <- performance(pred1,"tpr","fpr")
perf12 <- performance(pred1,"prec","rec")
pred2 <- prediction(res_rowmean$gkmSVM, res_rowmean$group)
perf21 <- performance(pred2,"tpr","fpr")
perf22 <- performance(pred2,"prec","rec")


pdf(file = "ROC2.pdf", width = 14, height = 7)
par(mfrow=c(1,2))

plot(perf11,xlab="1-Specificity",ylab="Sensitivity", lwd=2)
plot(perf21,colorize=F,lwd=2, col=2, add=T)
abline(coef = c(0,1), lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

plot(perf12,xlab="Recall",ylab="Precision", lwd=2, ylim=c(0,1))
plot(perf22,colorize=F,lwd=2, col=2, add=T)
abline(h=0.5, lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

dev.off()



# row max
precrec_obj3 <- evalmod(scores = res_rowmax$log_prob, labels = res_rowmax$group)
precrec_obj4 <- evalmod(scores = res_rowmax$gkmSVM, labels = res_rowmax$group)
precrec::auc(precrec_obj3)
precrec::auc(precrec_obj4)

pred1 <- prediction(res_rowmax$log_prob, res_rowmax$group)
perf11 <- performance(pred1,"tpr","fpr")
perf12 <- performance(pred1,"prec","rec")
pred2 <- prediction(res_rowmax$gkmSVM, res_rowmax$group)
perf21 <- performance(pred2,"tpr","fpr")
perf22 <- performance(pred2,"prec","rec")

pdf(file = "ROC2_max.pdf", width = 14, height = 7)
par(mfrow=c(1,2))

plot(perf11,xlab="1-Specificity",ylab="Sensitivity", lwd=2)
plot(perf21,colorize=F,lwd=2, col=2, add=T)
abline(coef = c(0,1), lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

plot(perf12,xlab="Recall",ylab="Precision", lwd=2, ylim=c(0,1))
plot(perf22,colorize=F,lwd=2, col=2, add=T)
abline(h=0.5, lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

dev.off()





# ------- may 21st -------
# motif freq ---- not 10-mer PWM score

# x=list_neg_1000[1:5,];lens=244
table_weight_prob = function(x, lens=350){
  # make Granges
  psn <- makeGRangesFromDataFrame(x, ignore.strand=T)
  # get seq from granges
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  # turn to characters
  seq_chr <- as.character(seq)
  
  # t=seq_chr[1]
  match_weight <- function(t){
    seq_vector <- mapply(substr, t, 1:(lens-9), 10:lens, USE.NAMES=FALSE)
    seq_vector_r <- chartr("ATGC","TACG",seq_vector)
    
    # weights
    m1 = d[match(seq_vector,de[,1]),2]
    m2 = d[match(seq_vector_r,d[,1]),2]
    m1[which(is.na(m1))]=m2[which(is.na(m1))]
    val_weight <- m1
    
    return(val_weight)
  }
  
  # t=seq_chr[1]
  match_prob <- function(t, chain){
    seq_vector <- mapply(substr, t, 1:(lens-14), 15:lens, USE.NAMES=FALSE)
    seq_vector_r <- chartr("ATGC","TACG",seq_vector)
    
    # log prod of probability (motif based)
    # s = seq_vector[1]
    seq_prob <- function(s){
      vec = unlist(strsplit(s,""))
      m1 = diag(pr1[vec,])
      m2 = diag(pr2[vec,])
      val_prob <- c(log(prod(m1)),log(prod(m2)))
      # return(val_prob)
      return(max(val_prob))
    }
    
    # log_prob <- do.call(rbind,mapply(seq_prob, seq_vector, SIMPLIFY = FALSE))
    # colnames(log_prob) <- c("pos_chain", "neg_chain")
    # 
    # if (chain == "pos"){
    #   return(log_prob[,1])
    # }else{
    #   return(log_prob[,2])
    # }
    
    log_prob <- do.call(c,mapply(seq_prob, seq_vector, SIMPLIFY = FALSE))
    return(log_prob)
    
  }
  
  xlist_w <- mapply(FUN=match_weight, seq_chr, SIMPLIFY = FALSE)
  # xlist_p <- mapply(FUN=match_prob, seq_chr, SIMPLIFY = FALSE, MoreArgs = list(chain="pos"))
  xlist_p <- mapply(FUN=match_prob, seq_chr, SIMPLIFY = FALSE, MoreArgs = list(chain="pos"))
  xlist_weight <- do.call(rbind,xlist_w)
  xlist_prob <- do.call(rbind,xlist_p)
  # xlist_weight <- do.call(rbind,lapply(xlist, function(x){x[[1]]}))
  # xlist_prob <- do.call(rbind,lapply(xlist, function(x){x[[2]]}))
  
  return(list(xlist_weight, xlist_prob))
}

# negSet_1000
a = proc.time()
res_neg <- table_weight_prob(x=list_neg_1000[1:1000,], lens = 244)
proc.time() - a

a = proc.time()
res_pos <- table_weight_prob(x=list_pos_1000[1:1000,], lens = 244)
proc.time() - a

save(res_neg, file="res_neg_CTCF.rda")
save(res_neg, file="res_pos_CTCF.rda")


# res_rowmean
res_rowmean <- cbind(rep(0:1, each=1000),
                     rbind(do.call(cbind, lapply(res_neg, function(x){rowMeans(x,na.rm = TRUE)})),
                           do.call(cbind, lapply(res_pos, function(x){rowMeans(x,na.rm = TRUE)}))))
colnames(res_rowmean) <- c("group", "gkmSVM", "log_prob")
res_rowmean <- as.data.frame(res_rowmean)

# res_max
res_rowmax <- cbind(rep(0:1, each=1000),
                    rbind(do.call(cbind, lapply(res_neg, function(x){apply(x,1,function(t){max(t, na.rm = TRUE)})})),
                          do.call(cbind, lapply(res_pos, function(x){apply(x,1,function(t){max(t, na.rm = TRUE)})}))))
colnames(res_rowmax) <- c("group", "gkmSVM", "log_prob")
res_rowmax <- as.data.frame(res_rowmax)



# row means
precrec_obj1 <- evalmod(scores = res_rowmax$log_prob, labels = res_rowmax$group)
precrec_obj2 <- evalmod(scores = res_rowmax$gkmSVM, labels = res_rowmax$group)
precrec_obj3 <- evalmod(scores = res_rowmean$log_prob, labels = res_rowmean$group)
precrec_obj4 <- evalmod(scores = res_rowmean$gkmSVM, labels = res_rowmean$group)
precrec::auc(precrec_obj1)
precrec::auc(precrec_obj2)
precrec::auc(precrec_obj3)
precrec::auc(precrec_obj4)
pred1 <- prediction(res_rowmean$log_prob, res_rowmean$group)
perf11 <- performance(pred1,"tpr","fpr")
perf12 <- performance(pred1,"prec","rec")
pred1 <- prediction(res_rowmax$log_prob, res_rowmax$group)
perf11 <- performance(pred1,"tpr","fpr")
perf12 <- performance(pred1,"prec","rec")

pred2 <- prediction(res_rowmean$gkmSVM, res_rowmean$group)
perf21 <- performance(pred2,"tpr","fpr")
perf22 <- performance(pred2,"prec","rec")
# pred2 <- prediction(res_rowmax$gkmSVM, res_rowmax$group)
# perf21 <- performance(pred2,"tpr","fpr")
# perf22 <- performance(pred2,"prec","rec")


pdf(file = "ROC2_motif_meanprob_meanSVM.pdf", width = 14, height = 7)
par(mfrow=c(1,2))

plot(perf11,xlab="1-Specificity",ylab="Sensitivity", lwd=2)
plot(perf21,colorize=F,lwd=2, col=2, add=T)
abline(coef = c(0,1), lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

plot(perf12,xlab="Recall",ylab="Precision", lwd=2, ylim=c(0,1))
plot(perf22,colorize=F,lwd=2, col=2, add=T)
abline(h=0.5, lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

dev.off()

















##########################
# test on K562 JUND
##########################
library(BSgenome.Hsapiens.UCSC.hg19)
library(gkmSVM)
library(seqinr)
library(precrec)
library(ggplot2)
library(gridExtra)
library(ROCR)

dir <- "H:/10-mers/20210521_test_CTCF_on_different_ChIPseq/"
setwd(dir)

# read in JUND from a new dataset ENCFF337DKJ.bed
bed <- as.data.frame(read.table(paste0(dir, "ENCFF337DKJ.bed"),
                                header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

JUND_new_bed <- bed[,1:3]
JUND_new_bed <- JUND_new_bed[order(JUND_new_bed[,1]),]
write.table(JUND_new_bed,
            file=paste0(dir, "ENCFF337DKJtest.bed"), quote=F, sep="\t", row.names=F, col.names=F)

# generate Nullseq
# inputBedFN = paste0("ENCFF337DKJtest.bed")
# genNullSeqs(inputBedFN)
# inputBedFN = JUND_new_bed
genNullSeqs("H:/10-mers/20210521_test_CTCF_on_different_ChIPseq/ENCFF337DKJtest.bed",
            genomeVersion = 'hg19')



set_prep <- function(file, number,lens=350){
  # read neg Set
  negSet <- read.fasta(file)
  neg_list0 <- strsplit(names(negSet),"_")
  # remove fix term
  len <- do.call(c, lapply(neg_list0, length))
  # form bed format
  neg_bed <- do.call(rbind, neg_list0[which(len == 5)])
  # generate list by chr
  neg_list <- as.data.frame(neg_bed[,1:3])
  neg_list[,4] <- as.numeric(neg_list[,3]) - as.numeric(neg_list[,2])+1
  colnames(neg_list) <- c("seqnames", "start", "end", "width")
  
  if (file == "negSet.fa"){
    output <- neg_list[1:number,]
  }else{
    bed1 <- bed[which(bed$V3-bed$V2==lens),]
    bed2 <- bed1[nrow(bed1):1,]
    # output <- neg_list[match(bed1[1:1000,2]+1,as.numeric(neg_list$start)),]
    output <- neg_list[match(bed2[1:1000,2]+1,as.numeric(neg_list$start)),]
  }
  
  return(output)
}

# neg Set 1000
list_neg_1000 <- set_prep("negSet.fa",1000)

# pos Set 1000
list_pos_1000 <- set_prep("posSet.fa",1000,lens=280)

save(list_neg_1000, file = "list_neg_1000.rda")
save(list_pos_1000, file = "list_pos_1000.rda")



# weights/pwm score
load(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/JUND/Part1/d.rda"))
de <- as.matrix(d)
dir.pwm <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/"
pwmfn <- paste0(dir.pwm, "JUND", ".ppm.txt")

# prepare pwm(+) and pwm(-)
pwm = t(as.matrix(read.table(pwmfn)))
pwm = (pwm + 0.001)/colSums(pwm + 0.001)

rownames(pwm) = c("A","C","G","T")
prob1 = pwm
rownames(prob1)=c("A","C","G","T")
pr1=(prob1)
# complimentary one
pr2=pr1[4:1,ncol(pr1):1]
rownames(pr2)=c("A","C","G","T")
rm(pwm,prob1,pwmfn)


# motif freq ---- not 10-mer PWM score

# x=list_neg_1000[1:5,];lens=244
table_weight_prob = function(x, lens=350){
  # make Granges
  psn <- makeGRangesFromDataFrame(x, ignore.strand=T)
  # get seq from granges
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  # turn to characters
  seq_chr <- as.character(seq)
  
  # t=seq_chr[1]
  match_weight <- function(t){
    seq_vector <- mapply(substr, t, 1:(lens-9), 10:lens, USE.NAMES=FALSE)
    seq_vector_r <- chartr("ATGC","TACG",seq_vector)
    
    # weights
    m1 = d[match(seq_vector,d[,1]),2]
    m2 = d[match(seq_vector_r,d[,1]),2]
    m1[which(is.na(m1))]=m2[which(is.na(m1))]
    val_weight <- m1
    
    return(val_weight)
  }
  
  # t=seq_chr[1]
  match_prob <- function(t, chain){
    seq_vector <- mapply(substr, t, 1:(lens-14), 15:lens, USE.NAMES=FALSE)
    seq_vector_r <- chartr("ATGC","TACG",seq_vector)
    
    # log prod of probability (motif based)
    # s = seq_vector[1]
    seq_prob <- function(s){
      vec = unlist(strsplit(s,""))
      m1 = diag(pr1[vec,])
      m2 = diag(pr2[vec,])
      val_prob <- c(log(prod(m1)),log(prod(m2)))
      # return(val_prob)
      return(max(val_prob))
    }
    
    # log_prob <- do.call(rbind,mapply(seq_prob, seq_vector, SIMPLIFY = FALSE))
    # colnames(log_prob) <- c("pos_chain", "neg_chain")
    # 
    # if (chain == "pos"){
    #   return(log_prob[,1])
    # }else{
    #   return(log_prob[,2])
    # }
    
    log_prob <- do.call(c,mapply(seq_prob, seq_vector, SIMPLIFY = FALSE))
    return(log_prob)
    
  }
  
  xlist_w <- mapply(FUN=match_weight, seq_chr, SIMPLIFY = FALSE)
  # xlist_p <- mapply(FUN=match_prob, seq_chr, SIMPLIFY = FALSE, MoreArgs = list(chain="pos"))
  xlist_p <- mapply(FUN=match_prob, seq_chr, SIMPLIFY = FALSE, MoreArgs = list(chain="pos"))
  xlist_weight <- do.call(rbind,xlist_w)
  xlist_prob <- do.call(rbind,xlist_p)
  # xlist_weight <- do.call(rbind,lapply(xlist, function(x){x[[1]]}))
  # xlist_prob <- do.call(rbind,lapply(xlist, function(x){x[[2]]}))
  
  return(list(xlist_weight, xlist_prob))
}

# negSet_1000
a = proc.time()
res_neg <- table_weight_prob(x=list_neg_1000[1:1000,], lens = 280)
proc.time() - a

a = proc.time()
res_pos <- table_weight_prob(x=list_pos_1000[1:1000,], lens = 280)
proc.time() - a

save(res_neg, file="res_neg_JUND.rda")
save(res_neg, file="res_pos_JUND.rda")


# res_rowmean
res_rowmean <- cbind(rep(0:1, each=1000),
                     rbind(do.call(cbind, lapply(res_neg, function(x){rowMeans(x,na.rm = TRUE)})),
                           do.call(cbind, lapply(res_pos, function(x){rowMeans(x,na.rm = TRUE)}))))
colnames(res_rowmean) <- c("group", "gkmSVM", "log_prob")
res_rowmean <- as.data.frame(res_rowmean)

# res_max
res_rowmax <- cbind(rep(0:1, each=1000),
                    rbind(do.call(cbind, lapply(res_neg, function(x){apply(x,1,function(t){max(t, na.rm = TRUE)})})),
                          do.call(cbind, lapply(res_pos, function(x){apply(x,1,function(t){max(t, na.rm = TRUE)})}))))
colnames(res_rowmax) <- c("group", "gkmSVM", "log_prob")
res_rowmax <- as.data.frame(res_rowmax)



# row means
precrec_obj1 <- evalmod(scores = res_rowmax$log_prob, labels = res_rowmax$group)
precrec_obj2 <- evalmod(scores = res_rowmax$gkmSVM, labels = res_rowmax$group)
precrec_obj3 <- evalmod(scores = res_rowmean$log_prob, labels = res_rowmean$group)
precrec_obj4 <- evalmod(scores = res_rowmean$gkmSVM, labels = res_rowmean$group)
precrec::auc(precrec_obj1)
precrec::auc(precrec_obj2)
precrec::auc(precrec_obj3)
precrec::auc(precrec_obj4)



# ---------------------
# max_prob / mean_SVM
# ---------------------
pred1 <- prediction(res_rowmax$log_prob, res_rowmax$group)
perf11 <- performance(pred1,"tpr","fpr")
perf12 <- performance(pred1,"prec","rec")

pred2 <- prediction(res_rowmean$gkmSVM, res_rowmean$group)
perf21 <- performance(pred2,"tpr","fpr")
perf22 <- performance(pred2,"prec","rec")
# pred2 <- prediction(res_rowmax$gkmSVM, res_rowmax$group)
# perf21 <- performance(pred2,"tpr","fpr")
# perf22 <- performance(pred2,"prec","rec")


pdf(file = "ROC2_JUND_maxprob_meanSVM.pdf", width = 14, height = 7)
par(mfrow=c(1,2))

plot(perf11,xlab="1-Specificity",ylab="Sensitivity", lwd=2)
plot(perf21,colorize=F,lwd=2, col=2, add=T)
abline(coef = c(0,1), lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

plot(perf12,xlab="Recall",ylab="Precision", lwd=2, ylim=c(0,1))
plot(perf22,colorize=F,lwd=2, col=2, add=T)
abline(h=0.5, lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

dev.off()



# ---------------------
# mean_prob / mean_SVM
# ---------------------
pred1 <- prediction(res_rowmean$log_prob, res_rowmean$group)
perf11 <- performance(pred1,"tpr","fpr")
perf12 <- performance(pred1,"prec","rec")

pred2 <- prediction(res_rowmean$gkmSVM, res_rowmean$group)
perf21 <- performance(pred2,"tpr","fpr")
perf22 <- performance(pred2,"prec","rec")


pdf(file = "ROC2_JUND_meanprob_meanSVM.pdf", width = 14, height = 7)
par(mfrow=c(1,2))

plot(perf11,xlab="1-Specificity",ylab="Sensitivity", lwd=2)
plot(perf21,colorize=F,lwd=2, col=2, add=T)
abline(coef = c(0,1), lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

plot(perf12,xlab="Recall",ylab="Precision", lwd=2, ylim=c(0,1))
plot(perf22,colorize=F,lwd=2, col=2, add=T)
abline(h=0.5, lty=3, col="grey")
legend("bottomright",bty="n",cex=1,
       legend = c("PWM", "gkmSVM"), 
       col=c("black","red"), lwd=c(2,2))

dev.off()













