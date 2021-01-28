options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)
arguments = matrix(unlist(strsplit(args,"=")),ncol=2,byrow = T)

for (args_i in 1:length(args)) {
  assign(arguments[args_i,1],as.numeric(arguments[args_i,2]))
}


library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)


motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", "JUND", 
          "MAX", "NANOG", "POU5F1", "RAD21", "REST", "RFX5", "SIX5",
          "SRF", "STAT1", "TCF12", "USF1", "USF2", "YY1")

# motif = "CTCF"
motif = motifs[ii]

digit = 3

# local path
# dir.pwm <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/"
# dir.save2 <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/final_dat/"
# dir.save3 <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Top_ranked_Dat/"

dir.pwm <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/PWM/"
dir.save2 <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/final_dat/"
dir.save3 <- "/scratch/bioinfo2/yjin85/Labtest/07062020_Motifs_over_SNP_deltaSVM/Top_ranked_Dat/"


# file path for pwm files
pwmfn <- paste0(dir.pwm, motif, ".ppm.txt")

# prepare pwm(+) and pwm(-) to calculate delta-PWM
pwm = t(as.matrix(read.table(pwmfn)))
# pwm = pwm + 0.01
# if(length(which(pwm==0))!=0){
#   pwm = pwm + 0.1
# }
pwm = (pwm + 0.001)/colSums(pwm + 0.001)

rownames(pwm) = c("A","C","G","T")
prob1 = pwm
rownames(prob1)=c("A","C","G","T")
pr1=log(prob1)
# complimentary one
pr2=pr1[4:1,ncol(pr1):1]
rownames(pr2)=c("A","C","G","T")
rm(pwm,prob1,pwmfn)

# motif length
motif.length = ncol(pr1)



# read in certain final_dat
load(paste0(dir.save2,"dat_",motif,".rda"))

# the dat are ranked by absolute deltaSVM values
dat.ranked <- dat[order(abs(dat$deltaSVM_19bp),decreasing = T),]

# the top 100 results
subdat <- dat.ranked[1:100,]

# get sequence from the table
ref.seq <- as.character(subdat$motif_ref)
mut.seq <- as.character(subdat$motif_mut)

# PWM score for top-ranked motif
PWM.ref <- mapply(function(x){sum(diag(pr1[unlist(strsplit(x,"")),1:motif.length]))},ref.seq)
PWM.mut <- mapply(function(x){sum(diag(pr1[unlist(strsplit(x,"")),1:motif.length]))},mut.seq)


dataset <- cbind(Chromosome = as.character(subdat$first.seqnames),
                 motif_start = subdat$second.X.start,
                 motif_end = subdat$second.X.end,
                 SNP_loc = subdat$first.start - subdat$second.X.start + 1,
                 motif_strand = as.character(subdat$second.X.strand),
                 motif_ref_seq = ref.seq,
                 motif_mutate_seq = mut.seq,
                 motif_ref_PWMscore = round(PWM.ref, digit),
                 motif_mutate_PWMscore = round(PWM.mut, digit),
                 deltaSVMsocre = round(subdat$deltaSVM_19bp,digit))


write.csv(dataset, file = paste0(dir.save3,"Table_",motif,".csv"),row.names = F)




