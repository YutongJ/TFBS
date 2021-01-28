options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)
arguments = matrix(unlist(strsplit(args,"=")),ncol=2,byrow = T)

for (args_i in 1:length(args)) {
  assign(arguments[args_i,1],as.numeric(arguments[args_i,2]))
}


library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

# # transform the pwm data to a easy-use version
# 
# for (i in c("GABPA","RAD21","SIX5")){
#   temp = read.table(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/Factorbook_",i,".txt"))
#   temp2 = base::rowSums(temp)
#   
#   write.table(signif(temp/temp2,6), paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/",i,".ppm.txt"),
#               row.names = F, col.names = F)
# }
# 
# i="EGR1"
# for (i in c("EGR1","JUND","MAX","REST","RFX5","SRF","USF1","USF2","YY1")){
#   temp = read.csv(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/JASPAR_",i,".csv"),header = F)
#   temp2 = base::colSums(temp)
#   
#   write.table(t(signif(temp/temp2,6)), paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/",i,".ppm.txt"),
#               row.names = F, col.names = F)
# }
# 
# rm(i,temp,temp2)



motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", "JUND", 
          "MAX", "NANOG", "POU5F1", "RAD21", "REST", "RFX5", "SIX5",
          "SRF", "STAT1", "TCF12", "USF1", "USF2", "YY1")


# BiocManager::install("rtracklayer")
# library(rtracklayer)

# write.table(tmp2,col.names = F,row.names = F,quote = F,
#             file = paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/TFBS.region/motif.",motif[i],".txt"))


# session <- browserSession()
# genome(session) <- "hg19"
# 
# query <- ucscTableQuery(session, "dbSnp153Composite", names = c("rs71375314"))
# 
# 
# tableName(query) <- "dbSnp153Common"
# test = getTable(query)




# motif = "CTCF"
motif = motifs[ii]

# # local paths
# dir.region <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/TFBS.region/"
# dir.d <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_mat/"
# dir.pwm <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/"
# dir.save1 <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/TF_over_SNPs/"
# dir.save2 <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/final_dat/"

# cluster paths
dir.region <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/TFBS.region/"
dir.d <- paste0("/scratch/bioinfo2/yjin85/TFBS/",motif,"/Part1/")
dir.pwm <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/PWM/"
dir.snp <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/SNPs_db153Common/"
dir.save1 <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/TF_over_SNPs/"
dir.save2 <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM//final_dat/"




# read in certain motif region on hg19
load(paste0(dir.region, "motif.",motif,".rdata"))
load(paste0(dir.d,"d.rda"))
de = as.matrix(d)

# file path for pwm files
pwmfn <- paste0(dir.pwm, motif, ".ppm.txt")

# prepare pwm(+) and pwm(-) to calculate delta-PWM
pwm = t(as.matrix(read.table(pwmfn)))
rownames(pwm) = c("A","C","G","T")
prob1 = pwm
rownames(prob1)=c("A","C","G","T")
pr1=log(prob1)
# complimentary one
pr2=pr1[4:1,ncol(pr1):1]
rownames(pr2)=c("A","C","G","T")
rm(pwm,prob1)

# motif length
motif.length = ncol(pr1)


# chr_num = 1
motif_overlap_chr <- function(chr_num){
  chr1 <- read.delim(paste0(dir.snp,"snp_chr",chr_num,".txt"),stringsAsFactors = F)
  chr1$alts <- substr(chr1$alts,1,1)
  
  # SNPs location: chr chromEnd
  chr1.tmp <- chr1[-1,c(1,3)]
  
  # prepare SNP data for "%over%" function with potential motif region
  chr1.snp.grange <- makeGRangesFromDataFrame(data.frame(seqnames=chr1.tmp[,1],
                                                         start=chr1.tmp[,2],
                                                         end=chr1.tmp[,2]), keep.extra.columns=T)
  # the indicator of which TF motif including SNPs
  # ind <- which( motif.TF %over% chr1.snp.grange)
  # ind2 <- which(chr1.snp.grange %within% motif.TF)
  # 1. two numbers are different because some motif regions are too closed with each other
  # examples: the motif are too closed to have overlapping region
  # chr1  1026065  1026079  15  - 0.8228106
  # chr1  1026068  1026082  15  - 0.7935567
  # 2. one motif can contain several SNPs
  
  
  # # ordered by chromStart
  # t1 <- tmp[ind,]
  # t1 <- t1[order(t1$start),]
  # t2 <- as.data.frame(chr1.snp.grange)[ind2,]
  
  
  
  # find all pairs between snps and motif containing certain snps
  pair1 <- findOverlapPairs(chr1.snp.grange,motif.TF)
  pair1 <- as.data.frame(pair1)
  snp1 <- chr1[match(pair1$first.end,chr1$chromEnd),5]
  snp2 <- chr1[match(pair1$first.end,chr1$chromEnd),6]
  pair1 <- cbind(pair1,snp1,snp2)
  return(list(pair1))
}

# test <- mapply(motif_overlap_chr, chr_num = 21:22)
# test <- motif_overlap_chr("21")

# motif_over_chrall <- do.call("rbind", mapply(motif_overlap_chr, chr_num = c(1:22,"X")))
# save(motif_over_chrall,file = paste0(dir.save1,"TF_chrs_",motif,".rda"))
# write.csv(motif_over_chrall, file = paste0(dir.save1,"TF_chrs_",motif,".csv"), row.names = F)


load(paste0(dir.save1,"TF_chrs_",motif,".rda"))
pairs <- motif_over_chrall

# get reference motif sequence in hg19 (with ref allele)
ref.seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,makeGRangesFromDataFrame(pairs[6:8], ignore.strand=T))
ref.seq <- as.character(ref.seq)

# find all corresponding snp
mut.snp <- pairs$snp2

loc <- pairs$first.end - pairs$second.X.start + 1

# verification
len <- length(which(substr(ref.seq,start = loc,stop = loc)!=pairs$snp1))
sink(file = paste0(dir.save1,"disparity_",motif,"_",len,".txt"))
which(substr(ref.seq,start = loc,stop = loc)!=pairs$snp1)
sink()

# get the mutate motif sequence in hg19
ref.mut <- paste0(substr(ref.seq,start = 1,stop = loc-1),mut.snp,substr(ref.seq,start = loc+1, stop = motif.length))




# sum(diag(pr1[unlist(strsplit(ref.seq[1],"")),1:motif.length]))
# sum(diag(pr1[unlist(strsplit(ref.mut[1],"")),1:motif.length]))
# 
# 
# test1[which(test1==-Inf)]=floor(min(test1[which(test1!=-Inf)]))



# prepare deltaSVM 19bps sequence
SNP_19bp <- data.frame(seqnames = pairs$first.seqnames,
                       start = pairs$first.end-9,
                       end = pairs$first.end+9)
# get 19bp ref sequence
snp_19bp_ref <- getSeq(BSgenome.Hsapiens.UCSC.hg19,makeGRangesFromDataFrame(SNP_19bp, ignore.strand=T))
snp_19bp_ref <- as.character(snp_19bp_ref)
# get 19bp alt sequence
snp_19bp_alt <- paste0(substr(snp_19bp_ref,start = 1,stop = 10-1),mut.snp,substr(snp_19bp_ref,start = 10+1, stop = 19))


# calculate gkmSVM from 19bps 
seq10_deltasum = function(x){
  xlist = mapply(substr, x, 1:10, 10:19, USE.NAMES=FALSE)
  xlistr= chartr("ATGC","TACG",xlist)
  
  m1 = d[match(xlist,de[,1]),2]
  m2 = d[match(xlistr,de[,1]),2]
  m1[which(is.na(m1))]=m2[which(is.na(m1))]
  sum = sum(m1)
  return(sum)
}
gkmSVM_19bp_ref = mapply(seq10_deltasum,snp_19bp_ref)
gkmSVM_19bp_alt = mapply(seq10_deltasum,snp_19bp_alt)
deltaSVM_19bp = gkmSVM_19bp_alt - gkmSVM_19bp_ref


# the overall table
dat <- data.frame(pairs[,c(1,2,7,8,10)],
                  motif_ref = ref.seq, motif_mut = ref.mut,
                  snp_19bp_ref,gkmSVM_19bp_ref,snp_19bp_alt,gkmSVM_19bp_alt,deltaSVM_19bp)

save(dat, file = paste0(dir.save2,"dat_",motif,".rda"))
write.csv(dat, file = paste0(dir.save2,"dat_",motif,".csv"),row.names = F)



