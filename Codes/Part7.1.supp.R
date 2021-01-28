options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

library("Biostrings")
wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
load(paste0(wd,"seq.ref.rda"))
load(paste0(wd,"seq.snp.rda"))
name = c("BCL11A","NANOG","JUN","POU5F1","STAT1","CTCF-K562","CTCF-H1")
span = c(11,11,13,11,11,15,15)

random.cluster <- function(i){
  TF = name[i]
    wd7 = paste0(wd,TF,"/Part7/")
    if(!dir.exists(wd7))  dir.create(wd7)
    load(paste0(wd,TF,"/Part1/d.rda"))
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
    
    gkmSVM.ref = mapply(seq10_deltasum,seq.ref)
    gkmSVM.snp = mapply(seq10_deltasum,seq.snp)
    save(gkmSVM.ref, file = paste0(wd7,"gkmSVM.ref.rda"))
    save(gkmSVM.snp, file = paste0(wd7,"gkmSVM.snp.rda"))
}

ii = as.numeric(args[1])
random.cluster(ii)