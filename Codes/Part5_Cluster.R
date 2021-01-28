##################################
###########  part 5  ############# summary of each position 
##################################


options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

p5 <- function(ii){
  name = c("CTCF","CTCF-H1","CTCF-K562")
  TF = name[ii]
  wd = paste0("/scratch/bioinfo2/yjin85/TFBS/")
  wd3.6 = paste0(wd,TF,"/Part3.6/")
  wd5 = paste0(wd,TF,"/Part5/")
  if(!dir.exists(wd5))  dir.create(wd5)
  
  summary.bps = function(loc1,num){
    top.loc = list()
    for(i in 1:20){
      fx = paste0(wd3.6,"/results.top",i,"/",loc1,"/top",i,loc1,seq(1:num),".alt3.rda")
      top1 = NULL
      for(j in 1:num){
        load(fx[j])
        out2 = as.vector(summary(as.vector(out)))
        top1 = rbind(top1,out2[1:6])
        # top1 = rbind(top1,as.vector(summary(as.vector(out))))
      }
      colnames(top1) = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
      top.loc[[i]] = top1
    }
    return(top.loc)
  }
  # loc1="motif";loc2="seqs";num=15
  top.pre.all =summary.bps(loc1 = "pre",num = 10)
  top.motif.all = summary.bps(loc1 = "motif",num = 15)
  top.aft.all = summary.bps(loc1 = "aft",num = 10)
  top.sum.all = list()
  for(i in 1:20){
    x = rbind(top.pre.all[[i]],top.motif.all[[i]],top.aft.all[[i]])
    top.sum.all[[i]]=x
  }
  save(top.pre.all,file = paste0(wd5,"top.pre.all.rda"))
  save(top.motif.all,file=paste0(wd5,"top.motif.all.rda"))
  save(top.aft.all,file = paste0(wd5,"top.aft.all.rda"))
  save(top.sum.all,file = paste0(wd5,"top.sum.all.rda"))
  
  wd3.6 = paste0(wd,TF,"/Part3.6/")
  save.summary1 = function(x,loc2){
    for(i in 1:20){
      dir.create(paste0(wd3.6,"/results.top",i,"/summary/"))
      sn.pre = paste0(wd3.6,"/results.top",i,"/summary/top",i,".",loc2,".rda")
      wn.pre = paste0(wd3.6,"/results.top",i,"/summary/top",i,".",loc2,".csv")
      # sn.pre = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011818/top",i,".",loc2,".rda")
      # wn.pre = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011818/top",i,".",loc2,".csv")
      top.loc = x[[i]]
      save(top.loc,file = sn.pre)
      top.loc = x[[i]]
      write.csv(top.loc,file = wn.pre,row.names = F)
    }
  }
  save.summary1(top.pre.all,loc2 = "pre")
  save.summary1(top.aft.all,loc2 = "aft")
  save.summary1(top.motif.all,loc2 = "motif")
  save.summary1(top.sum.all,loc2 = "sum")
  
  
  
  # get deltaSVM value matrix
  get.delta.bps = function(loc,num){
    top.loc = list()
    for(i in 1:20){
      fx = paste0(wd3.6,"/results.top",i,"/",loc,"/top",i,loc,seq(1:num),".alt3.rda")
      top1.loc = NULL
      for(j in 1:num){
        load(fx[j])
        top1.loc = cbind(top1.loc,out)
      }
      top.loc[[i]] = top1.loc
    }
    return(top.loc)
  }
  
  top.loc10pre = get.delta.bps(loc = "pre",num = 10)
  top.loc15motif = get.delta.bps(loc = "motif",num = 15)
  top.loc10aft = get.delta.bps(loc = "pre",num = 10)
  top.loc35 = list();for(i in 1:20){x = cbind(top.loc10pre[[i]],top.loc15motif[[i]],top.loc10aft[[i]]);top.loc35[[i]]=x}
  save(top.loc10pre,file = paste0(wd5,"/top.loc10pre.rda"))
  save(top.loc15motif,file = paste0(wd5,"/top.loc15motif.rda"))
  save(top.loc10aft,file = paste0(wd5,"/top.loc10aft.rda"))
  save(top.loc35,file = paste0(wd5,"/top.loc35.rda"))
}
ii = as.numeric(args[1])
p5(ii)