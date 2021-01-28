##################################
#########  part 3.4  ############# create ctcf datasets for cluster matching
##################################

# Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
# Input:
# x : randomseq from function "get.random()" / ctcfseq from function "get.ctcf()"
# len:  Length of motif
# fn: The lacation to store the result
# group: The motif type - random/ctcf

# Output:
# three prepared dataset: motif33, pre28, aft28 and their splitted sequences
options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

name = c("CTCF.old","REST")
span2=c(15,21)
p3.4<- function(i){
  TF = name[i]
  span = span2[i]
  wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
  # wd = paste0("D:/大�?�课程/emory/1-TFBS/结果/")
  wd3 = paste0(wd,TF,"/Part3/")
  wd3.4 = paste0(wd,TF,"/Part3.4/")
  if(!dir.exists(wd3.4))dir.create(wd3.4)
  # load(paste0(wd3,"motif.",TF,".rda"))
  
  
  get.TF <- function(pwmfn,fn){
    library(Biostrings)
    #### pwm downloaded from your folder
    pwm = t(as.matrix(read.table(pwmfn)))
    rownames(pwm) = c("A","C","G","T")
    
    # get reverse-comlementary seq
    pwm2=reverseComplement(pwm)
    
    #### get matched results
    library(BSgenome.Hsapiens.UCSC.hg19)
    motif.score=NULL
    chr.list=paste0("chr", c(1:22,"X"))
    for(chr in chr.list){
      message(chr)
      seq=Hsapiens[[chr]]
      hit1=matchPWM(pwm, unmasked(seq), with.score=T)
      if(length(hit1)>0){
        s1= mcols(hit1)$score
        s1= (s1-minScore(pwm))/(maxScore(pwm)-minScore(pwm))
        gr1=GRanges(seqnames=chr, ranges=hit1@ranges, strand="+", score=s1)
      }else{
        gr1=GRanges()
      }
      
      hit2=matchPWM(pwm2, unmasked(seq), with.score=T)
      if(length(hit2)>0){
        s2=hit2@elementMetadata$score
        s2= (s2-minScore(pwm2))/(maxScore(pwm2)-minScore(pwm2))
        gr2=GRanges(seqnames=chr, ranges=hit2@ranges, strand="-", score=s2)
      }else{
        gr2=GRanges()
      }
      
      gr=sort(c(gr1,gr2))
      motif.score=rbind(motif.score, as.data.frame(gr))
      motif.TF=makeGRangesFromDataFrame(motif.score, keep.extra.columns=T)
    }
    # rm(chr,gr,gr1,gr2,hit1,hit2,s1,s2,seq,motif.score,chr.list)
    rm(chr,gr,gr1,gr2,hit1,hit2,s1,s2,seq,chr.list)
    save(motif.score, file = paste0(fn,"motif.score.rda"))
    save(motif.TF, file = paste0(fn,"motif.",TF,".rda"))
  }
  get.TF(pwmfn = paste0(wd,TF,".ppm.txt"),fn=wd3)
  
  
  load(paste0(wd3,"motif.score.rda"))
  load(paste0(wd3,"motif.",TF,".rda"))
  
  get.datprep<-function(x, len, fn, group){
    dir.create(paste0(fn,"Seq.ref/"))
    dir.create(paste0(fn,"Seq.ref/list.motif",len,"/"))
    dir.create(paste0(fn,"Seq.ref/list.pre10/"))
    dir.create(paste0(fn,"Seq.ref/list.aft10/"))
    
    #get corresponding 33bps for 15-digit motif   33 = 15+9+9
    y = cbind(x[,1],x[,2]-9,x[,3]+9,x[,4:5])
    colnames(y) = colnames(x)[1:5]
    y[,4] = 9*2 + len
    y = makeGRangesFromDataFrame(y, keep.extra.columns=T)
    save(y, file = paste0(fn,"extend.motif",9*2+len,group,".rdata"))
    
    #get the corresponding 19bps for 15-digit motif  (33 in total)
    z = as.data.frame(y)
    for(i in 1:len){
      start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
      colnames(start.p) = colnames(z)
      psn = makeGRangesFromDataFrame(start.p, ignore.strand=T)
      seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
      seq_list = strsplit(as.character(seq),"")
      save(seq_list,file = paste0(fn,"Seq.ref/list.motif",len,"/","motif",i,".rda"))
    }
    rm(i,seq,psn,start.p,z)
    
    #get the corresponding 19bps for 10-digit pre motif (28 in total)
    pre28 = cbind(x[,1],x[,2]-19,x[,2]+8,x[,4:5])
    colnames(pre28) = colnames(x)[1:5]
    pre28[,4] = 28
    pre28 = makeGRangesFromDataFrame(pre28, keep.extra.columns=T)
    save(pre28, file = paste0(fn,"extend.pre28",group,".rdata"))
    
    #get the corresponding 19bps for 10-digit pre-motif
    z = as.data.frame(pre28)
    for(i in 1:10){
      start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
      colnames(start.p) = colnames(z)
      psn = makeGRangesFromDataFrame(start.p, ignore.strand=T)
      seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
      seq_list = strsplit(as.character(seq),"")
      save(seq_list,file = paste0(fn,"Seq.ref/list.pre10/pre",i,".rda"))
    }
    rm(i,seq,psn,start.p,z)
    
    #get the corresponding 19bps for 10-digit aft motif
    aft28 = cbind(x[,1],x[,3]-8,x[,3]+19,x[,4:5])
    colnames(aft28) = colnames(x)[1:5]
    aft28[,4] = 28
    aft28 = makeGRangesFromDataFrame(aft28, keep.extra.columns=T)
    save(aft28, file = paste0(fn,"extend.aft28",group,".rdata"))
    
    #get the corresponding 19bps for 10-digit aft-motif
    z = as.data.frame(aft28)
    for(i in 1:10){
      start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
      colnames(start.p) = colnames(z)
      psn = makeGRangesFromDataFrame(start.p, ignore.strand=T)
      seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
      seq_list = strsplit(as.character(seq),"")
      save(seq_list,file = paste0(fn,"Seq.ref/list.aft10/aft",i,".rda"))
    }
    rm(i,seq,psn,start.p,z)
    #####
  }
  get.datprep(x = motif.score, len = span, fn = wd3.4,group = TF)
  
  
  
  
  ##################################
  #########  part 3.5  ############# top sequence (the results are based on pwm matrix and hg19 scanning)
  ##################################
  
  # Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
  # Input:
  # fn: The lacation to store the result
  # ctcf: The ctcf results required from function "get.ctcf()"
  # len:  Length of motif
  # num: How many top sequence will be analysized in further steps (we adapt 20 here)
  
  # Output:
  # three prepared dataset: motif33, pre28, aft28 and their splitted sequences
  
  wd3.5 = paste0(wd,TF,"/Part3.5/")
  if(!dir.exists(wd3.5))dir.create(wd3.5)
  
  get.topseq = function(fn, ctcf, num=20,len=15){
    # get the corresponding 15bps seq in hg19
    seq15.hg19= getSeq(BSgenome.Hsapiens.UCSC.hg19,ctcf)
    
    seq15 = as.vector(seq15.hg19)
    cat("The number of unique motif is: ",length(unique(seq15)))  # 24402
    seq15 = sort(seq15)
    count.freq15 = cbind(unique(seq15),as.numeric(table(seq15)))
    count.freq15 = count.freq15[order(-as.numeric(table(seq15))),]
    save(count.freq15,file = paste0(fn,"count.freq",len,".rda"))
    
    # top.seq generating
    dir.create(paste0(fn,"toprow/"))
    toprow = paste0(fn,"toprow/top.row",seq(1:num),".rdata")
    for(i in 1:num){
      x = count.freq15[i,1]
      top.row = which(seq15 == x)
      save(top.row,file = toprow[i])
    }
    # # to check the amount of show-up time for corresponding complementary sequence (always zero)
    # for(i in 1:20){
    #   x = count.freq15[i,1]
    #   y = chartr("ATGC","TACG",x)
    #   m1 = length(which(seq15 == x))
    #   m2 = length(which(seq15 == y))
    #   m = cbind(m1,m2)
    #   print(m)
    # }
  }
  
  get.topseq(fn = wd3.5,motif.TF,num=20,len=span)
}

ii = as.numeric(args[1])
p3.4(ii)