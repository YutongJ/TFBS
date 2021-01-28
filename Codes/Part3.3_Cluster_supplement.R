##################################
#########  part 3.3  ############# random cluster
##################################
options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

library("Biostrings")
random.cluster <- function(ii){
  name = c("BCL11A","NANOG","JUN","POU5F1","STAT1","CTCF-K562","CTCF-H1")
  span2 = c(11,11,13,11,11,15,15)
  span1 = span2+20
  span = rep(0,length(name)+1)
  for(l in 1:length(name)){
    span[l+1] = sum(span1[1:l])
  }
  # i : number of TF
  for(i in 1:length(name)){
    if(ii <= span[i+1]) break
  }
  # j : 1-pre, 2-motif, 3-aft
  # k : number in file
  if(ii-span[i] <= 10){ 
    j = "pre"
    k = ii-span[i]
  }
  else if(ii-span[i] > 10 & ii-span[i] <= 10+span2[i]){ 
    j = "motif"
    k = ii-span[i]-10
  }
  else if(ii-span[i] > 10+span2[i]){ 
    j = "aft"  
    k = ii-span[i+1]+10
  }

  TF = name[i]
  if(j != "motif") len = 10
  else len = span2[i]
  wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
  wd1 = paste0(wd, TF, "/", "Part1/")
  wd3.2 = paste0(wd, TF, "/", "Part3.2/")
  wd3.3 = paste0(wd, TF, "/", "Part3.3/")
  if(!dir.exists(wd3.3))  dir.create(wd3.3)
  wd3.3 = paste0(wd3.3,j,"/")
  if(!dir.exists(wd3.3))  dir.create(wd3.3)
  fn = paste0(wd3.2,"Seq.ref/list.",j,len,"/",j,k,".rda")
  sn = paste0(wd3.3,"rand.",j,k,".alt3.rda")
  sn2 = paste0(wd3.3,"rand.",j,k,".alt3.csv")
  load(paste0(wd1,"d.rda"))
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
  
  load(fn)
  sub_seq = seq_list
  output = lapply(sub_seq,seq.delta)
  out = matrix(unlist(output),ncol = 3,byrow = T)
  save(out,file = sn)
  write.csv(out,file = sn2,row.names = F)
}

ii = as.numeric(args[1])
random.cluster(ii)