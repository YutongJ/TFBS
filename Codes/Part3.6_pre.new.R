##################################
#########  part 3.6  ############# ctcf cluster
##################################

options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

library("Biostrings")
ctcf.Cluster <- function(jj){
  name = c("CTCF","CTCF-H1","CTCF-K562")
  j = ceiling(jj/20)
  ii = jj%%20+1
  TF = name[j]
  wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
  wd3.6 = paste0(wd,TF,"/Part3.6/")
  if(!dir.exists(wd3.6)) dir.create(wd3.6)
  if(j == 1) load(paste0(wd,TF,"/Part1/d.rdata"))
  else load(paste0(wd,TF,"/Part1/d.rda"))
  load(paste0(wd,TF,"/Part3.5/toprow/top.row",ii,".rdata"))
  wd3.6 = paste0(wd3.6,"/results.top",ii,"/")
  if(!dir.exists(wd3.6)) dir.create(wd3.6)
  de = as.matrix(d)
  
  seq10_deltasum = function(x){
    xlist = c(paste(x[1:10 ], collapse = ''),paste(x[2:11 ], collapse = ''),paste(x[3:12 ], collapse = ''),paste(x[4:13 ], collapse = ''),paste(x[5:14 ], collapse = ''),
              paste(x[6:15 ], collapse = ''),paste(x[7:16 ], collapse = ''),paste(x[8:17 ], collapse = ''),paste(x[9:18 ], collapse = ''),paste(x[10:19 ], collapse = ''))
    xlistr= chartr("ATGC","TACG",xlist)
    
    m1 = d[match(xlist,de[,1]),2]
    m2 = d[match(xlistr,de[,1]),2]
    m1[which(is.na(m1))]=m2[which(is.na(m1))]
    sum = sum(m1)
    return(sum)
  }
  
  seq.delta = function(x){
    dict = c("A","T","G","C")
    temp = dict[which(dict!=x[10])]
    
    x1 = x; x2=x; x3=x
    x1[10]=temp[1]; x2[10]=temp[2]; x3[10]=temp[3]
    vec = NA
    vec[1] = seq10_deltasum(x1)-seq10_deltasum(x)
    vec[2] = seq10_deltasum(x2)-seq10_deltasum(x)
    vec[3] = seq10_deltasum(x3)-seq10_deltasum(x)
    return(vec)
  }
  
  wd3.4 = paste0(wd,TF,"/Part3.4/Seq.ref/")
  if(!dir.exists(paste0(wd3.6,"pre/"))) dir.create(paste0(wd3.6,"pre/"))
  fn = paste0(wd3.4,"list.pre10/pre",seq(1:10),".rda")
  sn = paste0(wd3.6,"pre/top",ii,"pre",seq(1:10),".alt3.rda")
  sn2= paste0(wd3.6,"pre/top",ii,"pre",seq(1:10),".alt3.csv")
  for(i in 1:10){
    print(i)
    load(fn[i])
    sub_seq = seq_list[top.row]
    output = lapply(sub_seq,seq.delta)
    out = matrix(unlist(output),ncol = 3,byrow = T)
    save(out,file = sn[i])
    write.csv(out,file = sn2[i],row.names = F)
  }
}

ii = as.numeric(args[1])
ctcf.Cluster(ii)