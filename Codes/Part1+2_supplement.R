##################################
###########  part 1  ############# Dateset preparation
##################################
options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

# variable name: TF, bedfn
# TF name: TF
# bed file name: bedfn
name = c("BCL11A","NANOG","JUN","POU5F1")
bedfile = c("ENCFF468MFY","ENCFF735WFO","ENCFF629BFI","ENCFF990CFV")


# input file path: wd
# output file path: outputwd
data.prep = function(i){
  TF = name[i]
  bedfn = bedfile[i]
  wd = "/scratch/bioinfo2/yjin85/TFBS/"
  outputwd = paste0(wd, TF, "/")
  if(!dir.exists(outputwd)) dir.create(outputwd)
  outputwd = paste0(outputwd, "Part1/")
  if(!dir.exists(outputwd)) dir.create(outputwd)
  
  # Package needed: seqinr, gkmSVM
  # Input:
  # wd: The path to store the posSet and negSet
  # dat: The name of original narrow peak ".bed" file downloaded from ENCODE
  # n: The number of posSet/negSet to be used in training steps
  
  # Output:
  # ".bed" file
  # posSet and negSet
  # 10000 randomly selected posSet and negSet (Notice: We use 10000 in our method, it can be other numbers)
  
  get.set<-function(fn,dat,n,ofn){
    setwd(fn)
    bed <- as.data.frame(read.table(paste0(dat,".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed1 = bed[,1:3]
    bed1 = bed1[order(bed1[,1]),]
    setwd(ofn)
    write.table(bed1, file=paste0(ofn, dat,"test.bed"), quote=F, sep="\t", row.names=F, col.names=F)
    # get original corresponding number of posSet and negSet
    library(gkmSVM)
    inputBedFN = paste0(dat,"test.bed")
    genNullSeqs(inputBedFN, genomeVersion = 'hg19')
    
    num = rep(NA,2)
    library(seqinr)
    test = read.fasta("posSet.fa")
    # pos
    int = sort(sample.int(length(test),min(length(test),n),replace = F))
    x = test[int]
    write.fasta(sequences = x, names = names(x), nbchar = 200, file.out = paste0("pos",min(length(test),n),".fa"))
    num[1] = min(length(test),n)
    # neg
    test = read.fasta("negSet.fa")
    int = sort(sample.int(length(test),min(length(test),n),replace = F))
    y = test[int]
    write.fasta(sequences = y, names = names(y), nbchar = 200, file.out = paste0("neg",min(length(test),n),".fa"))
    num[2] = min(length(test),n)
    return(num)
  }
  
  num = get.set(wd,bedfn,10000,outputwd)
  posneg = unlist(num)
  
  
  
  # Package needed: gkmSVM
  # Input:
  # posname: The name of positive set (.fa file)
  # nfn: The name of negative set (.fa file)
  # fn: The lacation to retrieve the training set and store the result
  # testfn: The test sequence file (all 10-mers included)
  
  # Output:
  # Kernel file
  # 2 GABPA_svmtrain files
  # output file with scores for each test sequence
  
  gkm.train<-function(fn,testfn,outname,ofn){
    setwd(fn)
    library(gkmSVM)
    # computes kernel
    posfn    = paste0(ofn,"pos",posneg[1],".fa")
    negfn    = paste0(ofn,"neg",posneg[2],".fa")
    kernelfn = paste0(ofn,TF,"kernel.txt")
    gkmsvm_kernel(posfn, negfn, kernelfn)
    svmfnprfx = paste0(ofn,TF,"_svmtrain")
    gkmsvm_train(kernelfn,posfn,negfn,svmfnprfx)
    
    outfn_nr10 = paste0(ofn,outname,".txt")
    gkmsvm_classify(testfn, svmfnprfx, outfn_nr10)
  }
  
  gkm.train(fn = wd,testfn = "nr10mers.fa.txt",outname = paste0("12878_",TF,"_nr10"),ofn = outputwd)
  
  
  
  # Package needed: None
  # Input:
  # vec: The splited 11-bps sequence
  # prob1: The log-transformed PWM for vec
  # prob2: The log-transformed PWM for complimentary vec
  
  # Output:
  # the largest log-transformed PWM value while matching target sequence to PWM
  
  ### 10-mers probability calculation (based on PWM)
  tenmers_prob<-function(vec,prob1=pr1,prob2=pr2){
    log_sum = rep(NA,(ncol(prob1)-9)*2)
    for(ii in 1:(ncol(prob1)-9)){
      num_s = ii
      num_e = ii+9
      m = diag(prob1[vec,num_s:num_e])
      log_sum[ii] = sum(m)
    }
    
    for(ii in 1:(ncol(prob1)-9)){
      num_s = ii
      num_e = ii+9
      m = diag(prob2[vec,num_s:num_e])
      log_sum[ii+(ncol(prob1)-9)] = sum(m)
    }  
    prob_max = max(log_sum)
    return(prob_max)
  }
  
  # fm: Frequency matrix of GABPA, downloaded from FACTORBOOK
  ### create the PWM from fm, pr1 and pr2
  
  pwmfn=paste0(TF,".ppm.txt")
  pwm = t(as.matrix(read.table(pwmfn)))
  rownames(pwm) = c("A","C","G","T")
  prob1 = pwm
  rownames(prob1)=c("A","C","G","T")
  pr1=log(prob1)
  # complimentary one
  pr2=pr1[4:1,ncol(pr1):1]
  rownames(pr2)=c("A","C","G","T")
  rm(pwm,prob1)
  
  # Package needed: None
  # Input:
  # outputfn: The path of output file containing scores for each test sequence
  # fn: The lacation to store the result
  
  # Output:
  # dataframe with 3 columns: sequences, SVM weights, PWM scores
  
  seq.score<-function(svmfn,ofn){
    setwd(ofn)
    ### read in output file --> (tenmer + weight)
    dataset = read.table(svmfn, sep="\t")
    data=strsplit(as.character(dataset[,1]),"")
    
    g = lapply(data,tenmers_prob)
    test1 = unlist(g)
    test1[which(test1==-Inf)]=floor(min(test1[which(test1!=-Inf)]))
    # save(test1,file = "test1.rdata")
    d = cbind(dataset,test1)
    save(d,file = paste0(ofn,"d.rda"))
    return(d)
  }
  
  d = seq.score(paste0("12878_",TF,"_nr10.txt"),outputwd)
  
  
  ### correlation between SVM weights and PWM score (for all 10-mers)
  cor1 <- cor(d$V2,d$test1)
  cor1 <- round(cor1, digits = 5)
  
  ### correlation between SVM weights and PWM score (for top 1000 10-mers)
  test1 = d[,3]
  num = test1[order(test1,decreasing = T)][1000]
  dd = subset(d,test1>=num)
  cor2 <- cor(dd$V2,dd$test1)
  cor2 <- round(cor2, digits = 5)
  
  ### make plot of SVM weights and PWM score
  # Package needed: ggplot2
  # Output:
  # All 10-mers TF.png
  # Top 1000 10-mers TF.png
  
  library(ggplot2)
  x = min(d$test1)+0.1*(max(d$test1)-min(d$test1))
  y = max(d$V2)-0.1*(max(d$V2)-min(d$V2))
  a <- ggplot(d, aes(x = test1, y = V2))+geom_point()
  a = a+ labs(title = paste0("All 10-mers  ", TF), y = "SVM weights", x = "PWM score")
  a = a+annotate("text", x = x, y = y, label = paste0("cor = ", cor1))+theme_bw()
  ggsave(paste0("All 10-mers ",TF,".png"), a, dpi = 300, width=8, height=6)
  
  x = min(dd$test1)+0.1*(max(dd$test1)-min(dd$test1))
  y = max(dd$V2)-0.1*(max(dd$V2)-min(dd$V2))
  b <- ggplot(dd, aes(x = test1, y = V2))+geom_point()
  b = b + labs(title = paste0("Top 1000 10-mers  ", TF), y = "SVM weights", x = "PWM score")
  b = b+annotate("text", x = x, y = y, label = paste0("cor = ", cor2))+theme_bw()
  ggsave(paste0("Top 1000 10-mers ", TF,".png"), b, dpi = 300, width=8, height=6)
  
  df = data.frame(TF,cor1,cor2)
  write.table(df,paste0(wd,"PWM-SVM cor.txt"), append = T, row.names = F,col.names = F)
  
  wd2 = paste0(wd, TF, "/", "Part2/")
  if(!dir.exists(wd2)) dir.create(wd2)
  
  # Package needed: None
  # Input:
  # dfn: The path of results from function "seq.score()" (.rda)
  # fn: The lacation to store the result
  
  # Output:
  # histogram
  # brief summary
  # return:  correlation between deltaPWM and deltaSVM for top 1000 sequence
  
  summary.var30<-function(dfn,fn){
    ### find top 1000 sequence
    load(dfn)
    de = as.matrix(d)
    # top 1000 sequence based on probability
    num = order(d[,3], decreasing = T)[1:1000]
    sub_de = d[num,]
    save(sub_de,file = paste0(fn,"sub_de_prob.rdata"))
    
    # top 1000 sequence based on gkmSVM value
    num = order(d[,2],decreasing = T)[1:1000]
    sub_de_gkm = de[num,]
    save(sub_de_gkm,file = paste0(fn,"sub_de_gkm.rdata"))
    
    sub_de = sub_de_gkm
    seq_list = strsplit(as.character(sub_de[,1]),"")
    
    # create single variant sequence
    seq_ex = function(x){
      data_seq = matrix(rep(x,30),ncol = 10,byrow = T)
      dict = c("A","T","G","C")
      for (i in 1:10) {
        temp <- dict[which(dict!=x[i])]
        for (j in 1:3) {
          data_seq[(i-1)*3+j,i] = temp[j]
        }
      }
      vec=rep(NA,30)
      for(ii in 1:30){
        vec[ii] = paste(data_seq[ii,], collapse = '')
      }
      return(vec)
    }
    m = lapply(seq_list,seq_ex)
    
    # create single variant sequence matrix (30 variant in one row)------ n
    n = matrix(unlist(m),ncol = 30,byrow = T)
    rm(num,seq_ex,m,seq_list)
    
    # create single variant reverse comlementary sequence------ n_c
    n_c1 = chartr("ATGC","TACG",n)
    n_c = n_c1
    for(i in 1:1000){
      for(j in 1:30){
        n_c[i,j]=paste(rev(strsplit(n_c1[i,j],"")[[1]]),collapse = '')
      }
    }
    remove(n_c1,i,j)
    
    # sub_de refers to top 1000 tenmers ranked by gkmSVM weights
    test1 = as.vector(t(n))
    test2 = as.vector(t(n_c))
    m1 = match(test1,d[,1])
    m2 = match(test2,d[,1])
    mm = m1
    mm[which(is.na(mm))]=m2[which(is.na(mm))]
    rm(test1,test2,m1,m2)
    
    # matched sequence for changed sequence
    x.seq = matrix(d[mm,1],ncol = 30,byrow = T)
    # original gkmSVM value for changed sequence
    dist.ori = matrix(d[mm,2],ncol = 30,byrow = T)
    # difference between changed sequence and original seq
    dist = matrix(NA,ncol = 30,nrow = 1000)
    for (i in 1:1000){
      for(j in 1:30){
        dist[i,j]=dist.ori[i,j]-as.numeric(sub_de_gkm[i,2])
      }
    }
    rm(i,j)
    
    # original probability for changed sequence based on probability
    dist.p.ori = matrix(d[mm,3],ncol = 30,byrow = T)
    # difference probability between changed sequence and original seq
    dist.p = matrix(NA,ncol = 30,nrow = 1000)
    for (i in 1:1000){
      for(j in 1:30){
        dist.p[i,j]=exp(dist.p.ori[i,j])-exp(as.numeric(sub_de_gkm[i,3]))
      }
    }
    rm(i,j)
    
    # correlation between deltaPWM and deltaSVM for top 1000 sequence
    g = as.vector(t(dist))
    g.p = as.vector(t(dist.p))
    cr=round(cor(g,g.p), digits = 5) 
    
    # plot for correlation between deltaPWM and deltaSVM for top 1000 sequence 
    library(ggplot2)
    xmin = min(g)+0.1*(max(g)-min(g))
    ymax = max(g.p)-0.1*(max(g.p)-min(g.p))
    d = cbind(data.frame(g),data.frame(g.p))
    b <- ggplot(d, aes(x = g, y = g.p))+geom_point()+theme_bw()
    b = b + labs(title = paste0("deltaPWM and deltaSVM for top 1000 sequences ", TF), y = "deltaPWM", x = "deltaSVM")
    b = b+annotate("text", x = xmin, y = ymax, label = paste0("cor = ", cr))
    ggsave(paste0(fn,"Delta Top 1000 10-mers ",TF,".png"), b, width=8, height=6)
    
    # histogram for all 30000 changed sequence
    g = unique(as.numeric(as.vector(t(dist))))
    jpeg(file=paste0(fn,"Distribution of single variant.jpeg"))
    hist(g,main = paste0("distribution of ",length(g)," single variant ", TF),xlab = "difference")
    dev.off()
    
    # summary for all 30000 changed sequence
    sink(file = paste0(fn,"summary of ",length(g)," single variant"))
    summary(g)
    sink()
    return(cr)
  }
  
  cr = summary.var30(dfn=paste0(wd1,"d.rda"),fn=wd2)
  df = data.frame(TF,cr)
  write.table(df,paste0(wd,"delta cor.txt"), append = T, row.names = F,col.names = F)
}

ii = as.numeric(args[1])
data.prep(ii)