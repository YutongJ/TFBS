# load packages
library(dplyr)
library(ggplot2)
library(gridExtra)
# library(RVenn)
library(extrafont)

# motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", 
#           "JUND", "MAX", "NANOG", "POU5F1", "RAD21", 
#           "REST", "RFX5", "SIX5", "SRF", "STAT1", 
#           "TCF12", "USF1", "USF2", "YY1")

#--------New: Remove REST
motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", 
          "JUND", "MAX", "NANOG", "POU5F1", "RAD21", 
          "RFX5", "SIX5", "SRF", "STAT1", 
          "TCF12", "USF1", "USF2", "YY1")

# for (j in 1:19) {
#   dir.create(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/",motifs[j]))
# }
# 
# for (j in 1:19) {
#   dir.create(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/",motifs[j],"/Part1"))
# }

# motif = motifs[2]
one_plot <- function(motif,num){
  load(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/",motif,"/Part1/d.rda"))
  
  ### correlation between SVM weights and PWM score (for top 1000 10-mers)
  dd <- d[order(d[,3],decreasing = T),]
  dd <- dd[1:1000,]
  
  cor2 <- round(cor(dd$V2,dd$test1), digits = 5)
  
  ### make plot of SVM weights and PWM score
  # Package needed: ggplot2
  # Output:
  # All 10-mers TF.png
  # Top 1000 10-mers TF.png
  
  x = max(dd$test1)-0.02*(max(d$test1)-min(d$test1))
  y = min(dd$V2)
  b <- ggplot(dd, aes(x = test1, y = V2))+geom_point()
  b = b + labs(title = paste0(letters[num],". ",motif), y = "SVM weights", x = "PWM score")
  b = b+annotate("text", x = x, y = y, label = paste0("cor = ", cor2))+theme_bw()
  return(b)
}

one_plot(motifs[4],num = 4)

panels = list()
for (j in seq(motifs)) {
  panels[[j]] <- one_plot(motifs[j],num=j)
}

temp = motifs[c(1,17,13,1)]
#--------New: Remove REST
temp = motifs[c(1,17,13,1)]
panels = list()
for (j in seq(temp)) {
  panels[[j]] <- one_plot(temp[j],num=j)
}

pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot11.pdf",width = 8)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],nrow = 2)
dev.off()


pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot1.pdf",width = 8)
grid.arrange(panels[[2]],panels[[17]],panels[[13]],panels[[1]],nrow = 2)
dev.off()


pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot22.pdf",width = 20,height = 14)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],panels[[5]],
             panels[[6]],panels[[7]],panels[[8]],panels[[9]],panels[[10]],
             panels[[11]],panels[[12]],panels[[13]],panels[[14]],panels[[15]],
             panels[[16]],panels[[17]],panels[[18]],panels[[19]],nrow = 4)
dev.off()

#--------New: Remove REST
pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Supplementary Figure1.pdf",width = 20,height = 14)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],panels[[5]],
             panels[[6]],panels[[7]],panels[[8]],panels[[9]],panels[[10]],
             panels[[11]],panels[[12]],panels[[13]],panels[[14]],panels[[15]],
             panels[[16]],panels[[17]],panels[[18]],nrow = 4)
dev.off()






# prepare excel
require(openxlsx)

tables = list()
for (i in seq(motifs)) {
  tmp <- read.csv(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Top_ranked_Dat/","Table_",motifs[i],".csv"),stringsAsFactors = F)
  tmp <- tmp[,-5]
  tmp <- cbind(tmp[,c(1,2,3,4,5,7,6,8)],tmp[,8]-tmp[,7],tmp[,9])
  colnames(tmp) <- c("Chromosome","motif_start","motif_end","SNP_loc","motif_ref_seq",       
                     "motif_ref_PWMscore","motif_alt_seq","motif_alt_PWMscore","motif_deltaPWM_score","deltaSVMsocre")
  tables[[i]] <- tmp
}

names(tables) <- motifs
write.xlsx(tables, file = paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/","Top100deltaSVM_of_SNP_in_Motifs.xlsx"),
           colWidths = rep("auto",length(motifs)))
#--------New: Remove REST
write.xlsx(tables, file = paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/","Supplementary Table1.xlsx"),
           colWidths = rep("auto",length(motifs)))





# prepare barplots
rm(list=ls())

motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", "JUND", 
          "MAX", "NANOG", "POU5F1", "RAD21", "REST", "RFX5", "SIX5",
          "SRF", "STAT1", "TCF12", "USF1", "USF2", "YY1")

#--------New: Remove REST
motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", "JUND", 
          "MAX", "NANOG", "POU5F1", "RAD21", "RFX5", "SIX5",
          "SRF", "STAT1", "TCF12", "USF1", "USF2", "YY1")

motif.lengths = c()
motif = motifs[2]
bar_plot <- function(motif){
  load(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/",motif,"/Part1/d.rda"))
  
  x = d[order(d[,2],decreasing = T),]
  x = as.matrix(x[1:20,])

  # prepare dataframe for plotting
  a <- rep(as.character(x[,1]),each=2)
  b <- rep(c("Weight","10-mer PWM score"),20)
  c <- rep(0,40)
  c[(2*(1:20)-1)] <- as.numeric(x[,2])
  c[(2*(1:20))] <- as.numeric(x[,3])
  plot_seq <- data.frame("Sequence"=a,
                         "Scores"=b,
                         "Freq"=c,stringsAsFactors = F)
  # plots
  p <- 
    plot_seq %>% 
    ggplot(aes(x = Sequence, y = Freq, group = Scores, fill = Scores)) +
    geom_bar(stat = "identity", width = 0.75) +
    coord_flip() +
    scale_x_discrete(limits = rev(unique(plot_seq[,1]))) +
    # another trick!
    # scale_y_continuous(breaks = seq(-15, 15, 3),
    #                    labels = abs(seq(-15,15,3))) +
    labs(x = "10-mer", y = "Values", title = motif) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill =  "grey90")) +
    # reverse the order of items in legend
    # guides(fill = guide_legend(reverse = TRUE)) +
    # change the default colors of bars
    scale_fill_manual(values=c("darkblue", "black"),
                      name="",
                      breaks=c("10-mer PWM score","Weight"),
                      labels=c("log-transformed 10-mer PWM score","Weight")) 
  return(p)
}

bar_plot(motifs[2])


panels = list()
for (j in seq(motifs)) {
  panels[[j]] <- bar_plot1(motifs[j])
}


pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot3.pdf",width = 8)
grid.arrange(panels[[2]],panels[[17]],panels[[13]],panels[[1]],nrow = 2)
dev.off()


pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot4.pdf",width = 20,height = 14)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],panels[[5]],
             panels[[6]],panels[[7]],panels[[8]],panels[[9]],panels[[10]],
             panels[[11]],panels[[12]],panels[[13]],panels[[14]],panels[[15]],
             panels[[16]],panels[[17]],panels[[18]],panels[[19]],nrow = 4)
dev.off()

# par(mfrow=c(1,1))
panels.old = panels







# prepare barplots
rm(list=ls())

motifs= c("BCL11A", "CTCF", "EGR1", "GABPA", "JUN", "JUND", 
          "MAX", "NANOG", "POU5F1", "RAD21", "REST", "RFX5", "SIX5",
          "SRF", "STAT1", "TCF12", "USF1", "USF2", "YY1")

motif.lengths = c(11, 15, 14, 11, 13, 11, 10, 11, 11, 15, 21, 15, 15, 18, 11, 11, 11, 16, 12)

# motif = motifs[1]
motif.length = motif.lengths[1]
dir.pwm <- "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/PWM/"


### 10-mers probability calculation (based on PWM)
# tenmer = as.character(y[1,1]);prob1=pr1;prob2=pr2
tenmers_prob<-function(tenmer,prob1=pr1,prob2=pr2){
  vec = unlist(strsplit(tenmer,""))
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



# scaled version
bar_plot1 <- function(motif,motif.length,num){
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
  # rm(pwm,prob1,pwmfn)
  
  load(paste0("H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/d_matrix/",motif,"/Part1/d.rda"))
  
  x = d[order(d[,2],decreasing = T),]
  # x = as.matrix(x[1:20,])
  z=x[1:20,]
  y=x[1:20,]
  y[,3] = mapply(tenmers_prob,as.character(y[,1]), MoreArgs = list(prob1=pr1,prob2=pr2))
  y[,3] = y[,3]/max(abs(y[,3]))
  y[,2] = y[,2]/max(abs(y[,2]))
  

  # prepare dataframe for plotting
  a <- rep(as.character(y[,1]),each=2)
  b <- rep(c("Weight","10-mer PWM score"),20)
  c <- rep(0,40)
  c[(2*(1:20)-1)] <- y[,2]
  c[(2*(1:20))] <- (y[,3])
  plot_seq <- data.frame("Sequence"=a,
                         "Scores"=b,
                         "Freq"=c,stringsAsFactors = F)
  # plots
  p <- 
    plot_seq %>% 
    ggplot(aes(x = Sequence, y = Freq, group = Scores, fill = Scores)) +
    geom_bar(stat = "identity", width = 0.75) +
    coord_flip() +
    scale_x_discrete(limits = rev(unique(plot_seq[,1]))) +
    # another trick!
    scale_y_continuous(breaks = seq(-1, 1, 1),
                       labels = c(min(round(z[,3],1)),0,max(round(z[,2],1)))) +
    labs(x = "10-mer", y = "Values", title = paste0(letters[num],". ",motif)) +
    theme(legend.position = "bottom",text=element_text(family="mono"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill =  "grey90")) +
    # reverse the order of items in legend
    # guides(fill = guide_legend(reverse = TRUE)) +
    # change the default colors of bars
    scale_fill_manual(values=c("darkblue", "black"),
                      name="",
                      breaks=c("10-mer PWM score","Weight"),
                      labels=c("log PWM score","Weight")) 
  return(p)
}

# bar_plot1(motifs[1],motif.lengths[1],num=1)

# panels = list()
# for (j in seq(motifs)) {
#   panels[[j]] <- bar_plot1(motifs[j],motif.lengths[j])
# }


# pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot32.pdf",width = 6)
# grid.arrange(panels[[2]],panels[[17]],panels[[13]],panels[[1]],nrow = 2)
# dev.off()



temp = motifs[c(1,17,13,1)]
temp2 = motif.lengths[c(1,17,13,1)]
panels = list()
for (j in seq(temp)) {
  panels[[j]] <- bar_plot1(temp[j],temp2[j],num=j)
}

pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot33.pdf",width = 8)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],nrow = 2)
dev.off()

panels = list()
for (j in seq(motifs)) {
  panels[[j]] <- bar_plot1(motifs[j],motif.lengths[j],num=j)
}

pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Plot42.pdf",width = 15,height = 14)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],panels[[5]],
             panels[[6]],panels[[7]],panels[[8]],panels[[9]],panels[[10]],
             panels[[11]],panels[[12]],panels[[13]],panels[[14]],panels[[15]],
             panels[[16]],panels[[17]],panels[[18]],panels[[19]],nrow = 4)
dev.off()
#--------New: Remove REST
pdf(file = "H:/10-mers/07062020_Motifs_over_SNP_deltaSVM/Paper_Supp/Supplementary Figure2.pdf",width = 15,height = 14)
grid.arrange(panels[[1]],panels[[2]],panels[[3]],panels[[4]],panels[[5]],
             panels[[6]],panels[[7]],panels[[8]],panels[[9]],panels[[10]],
             panels[[11]],panels[[12]],panels[[13]],panels[[14]],panels[[15]],
             panels[[16]],panels[[17]],panels[[18]],nrow = 4)
dev.off()



library(purrr)
set.seed(42)
toy = map(sample(5:25, replace = TRUE, size = 10),
          function(x) sample(letters, size = x))

toy <- Venn(toy)
ggvenn(toy,slice = c(2,3))





