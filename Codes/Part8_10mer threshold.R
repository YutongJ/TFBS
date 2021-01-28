# get the threshold for 10-mers

motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", 
          "JUND", "MAX", "NANOG", "POU5F1", "RAD21", 
          "REST", "RFX5", "SIX5", "SRF", "STAT1", 
          "TCF12", "USF1", "USF2", "YY1")

motif = motifs[2]

# the size of 
sizes <- 100000

# ten-mer threshold
load(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/",motif,"/Part1/d.rda"))
# de = as.matrix(d)

# test  = d
d[,1] = as.character(d[,1])
# 

# generate the random collection of 10-mers
d_rand <- d[sample(seq(nrow(d)), size = sizes, replace = F),]

# random positions within 10-mers
pos <- sample(1:10,size = sizes,replace = T)

# the reference 10-mer sequences
seq.ref <- as.character(d_rand[,1])
# the SNP in reference 10-mer
d_rand_ref <- substr(seq.ref,pos,pos)

# function to change targeted SNP within the 10-mer sequence
# x = d_rand_ref[1]
mut_allele <- function(x){
  temp <- c("A","C","G","T")
  temp2<- temp[-which(temp==x)]
  xx <- temp2[sample(1:3,1)]
  return(xx)
}

# get the SNP in the corresponding alternative 10-mers
d_rand_mut <- mapply(mut_allele,d_rand_ref)

# the alternative 10-mer sequences
seq.mut <- paste0(substr(seq.ref,start = 1,stop = pos-1),d_rand_mut,substr(seq.ref,start = pos+1, stop = 10))


# match corresponding weights for 10-mer
seq10_weight = function(x){
  xlistr= chartr("ATGC","TACG",x)
  
  m1 = d[match(x,d[,1]),2]
  m2 = d[match(xlistr,d[,1]),2]
  m1[which(is.na(m1))]=m2[which(is.na(m1))]
  return(m1)
}
# seq10_weight(seq.ref[1])



# get the weights for reference 10-mer
weight.ref <- d[match(seq.ref, d[,1]), 2]

seq.mut2 <- chartr("ATGC","TACG",seq.mut)

m1 <- d[match(seq.mut, d[,1]), 2]
m2 <- d[match(seq.mut2, d[,1]), 2]
m1[which(is.na(m1))]=m2[which(is.na(m1))]

# get the weights for alternative 10-mer
weight.mut <- m1
rm(m1,m2,seq.mut2)


# the weight difference of 10000 10-mers should be:
weight_diff <- weight.mut - weight.ref


# the empirical quantiles
quantile(weight_diff,probs = c(0.025,0.975))




# you can use sink() here to get the results saved in a txt file
# here is some example code, please revise according to your settings
sink(file = "Your path to save the results", append = T) # if you set append=T and use the same file name, it will append in the same file
cat(c(motif,quantile(weight_diff,probs = c(0.025,0.975))))
sink()






