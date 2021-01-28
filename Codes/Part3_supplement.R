##################################
#########  part 3.1  ############# prepare datasets for 15-bp ctcf and randomly selected sequence
##################################
options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

name = c("BCL11A","NANOG","JUN","POU5F1","STAT1","CTCF-K562","CTCF-H1")
span2 = c(11,11,13,11,11,15,15)
P3 <- function(ii){
  TF = name[ii]
  span = span2[ii]
  wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
  wd1 = paste0(wd, TF, "/", "Part1/")
  wd3 = paste0(wd, TF, "/", "Part3/")
  dir.create(wd3)
  
  # Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
  # Input:
  # pwmfn: The path of pwm (PWM in "A","C","G","T" row order)
  # fn: The lacation to store the result
  
  # Output:
  # all matched TF location on Genome
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
    rm(chr,gr,gr1,gr2,hit1,hit2,s1,s2,seq,motif.score,chr.list)
    save(motif.TF, file = paste0(fn,"motif.",TF,".rdata"))
  }
  get.TF(pwmfn = paste0(wd,TF,".ppm.txt"),fn=wd3)

# Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19, stringr
# Input:
# number: The number of regions to sample
# span: The nubmer of bp to add to start
# fn: The lacation to store the result
# ctcffn: The path of ctcf results required from function "get.ctcf()"

# Output:
# random selected sequence set.

  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(stringr)

  setwd(wd3)
  
  get.random<-function(tffn,fn,number,span){
    #chromosomes of interest
    my_chr <- c(1:22,'X')
    my_chr <- gsub(pattern="^", replacement='chr', my_chr)
    
    #initialise list to store chromosome sizes
    my_chr_size <- list()
    for (i in my_chr){
      my_chr_size[[i]] <- length(BSgenome.Hsapiens.UCSC.hg19[[i]])
    }
    
    #checkout my_chr_size
    head(my_chr_size,2)
    save(my_chr_size,file = paste0(fn,"my_chr_size.rda"))
    
    #initialise some vectors for storing random coordinates
    my_random_start  <- vector()
    my_random_end    <- vector()
    my_random_chr    <- vector()
    my_random_strand <- vector()
    
    #loop through number of regions
    for(i in 1:number){
      my_random_chr[i] <- sample(x=my_chr,size=1)
      my_random_strand[i] <- sample(x=c('-','+'),size=1)
      my_max <- my_chr_size[[my_random_chr[i]]]-span
      my_random_start[i] <- floor(runif(n=1, min=20, max=my_max-19))
      my_random_end[i] <- my_random_start[i] + span
    }
    my_random_width <- rep(span+1,number)
    my_random_seq <- data.frame(seqnames=my_random_chr,
                                start=my_random_start,
                                end=my_random_end,
                                width = my_random_width,
                                strand=my_random_strand)
    head(my_random_seq)
    rm(i,my_chr,my_random_chr,my_random_start,my_random_end,my_random_strand,number,span,my_max,my_random_width)
    random_seq = my_random_seq[order(my_random_seq$seqnames,my_random_seq$start),]
    
    load(tffn)
    random_seq1 = makeGRangesFromDataFrame(random_seq, keep.extra.columns=T)
    x = random_seq1 %over% motif.TF
    randomseq = random_seq[-which(x),]
    save(randomseq,file = paste0(fn,"randomseq.rda"))
  }
  get.random(tffn=paste0(wd3,"motif.",TF,".rdata"),fn = wd3,number=10000,span = span-1)




##################################
#########  part 3.2  ############# create random datasets for cluster matching
##################################

# Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
# Input:
# x : randomseq from function "get.random()" / ctcfseq from function "get.ctcf()"
# len:  Length of motif
# fn: The lacation to store the result
# group: The motif type - random/ctcf

# Output:
# three prepared dataset: motif33, pre28, aft28 and their splitted sequences

  wd3.2 = paste0(wd, TF, "/", "Part3.2/")
  dir.create(wd3.2)
  setwd(wd1)
  span = ncol(t(read.table(paste0(wd,TF,".ppm.txt"))))
  # pos.set <- Sys.glob("pos*.fa")
  # pos.num = str_extract_all(pos.set[1],"\\d")
  # pos.num = unlist(pos.num[[1]])
  # pos.num = as.numeric(paste(pos.num[1:length(pos.num)], collapse = ""))
  load(paste0(wd3,"randomseq.rda"))
  
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
  get.datprep(x = randomseq, len = span, fn = wd3.2,group = "random")
}

ii = as.numeric(args[1])
P3(ii)