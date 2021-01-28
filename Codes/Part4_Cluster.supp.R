##################################
###########  part 4  ############# threshold
##################################
options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
name = c("BCL11A","NANOG","JUN","POU5F1","STAT1","CTCF-K562","CTCF-H1")
span = c(11,11,13,11,11,15,15)

# Package needed: None
# Input:
# fn: The path of first cell line pre/motif/aft sequences
# Output:
# return:  threshold for prediction based on three cell line

#####  random seq
get.cre<-function(i){
  TF = name[i]
  len = span[i]
  wd7 = paste0(wd, TF, "/Part7/")
  if(!dir.exists(wd7)) dir.create(wd7)
  sum.bps = function(fn,TF,loc,num){
    fx = paste0(fn,loc,"/rand.",loc,seq(1:num),".alt3.rda")
    top1 = NULL
    for(j in 1:num){
      load(fx[j])
      top1 = cbind(top1,as.vector(out))
    }
    return(top1)
  }
  
  cell.all = function(x){
    r1 = sum.bps(fn = paste0(wd,TF,"/Part3.3/"),loc = "pre",num = 10)
    r2 = sum.bps(fn = paste0(wd,TF,"/Part3.3/"),loc = "motif",num = len)
    r3 = sum.bps(fn = paste0(wd,TF,"/Part3.3/"),loc = "aft",num = 10)
    rand.all = cbind(r1,r2,r3)
    return(rand.all)
  }
  
  rand.gm  = cell.all()
  # cre = quantile(as.vector(rand.gm),probs=0.05,na.rm = T)
  save(rand.gm,file = paste0(wd7,"rand.gm.rda"))
}

ii = as.numeric(args[1])
get.cre(ii)