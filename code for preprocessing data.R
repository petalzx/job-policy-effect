setwd("/Users/cq62/Library/CloudStorage/Dropbox/Research/KitagawaLeeQiu_personal/ML connections/empirical app-JTPA-new/no_cost")
rm(list=ls())

jtpa.cq <- read.csv("mydata.csv",header = TRUE)[,-(1:2)]
#jtpa.cq <-jtpa.cq[-which(jtpa.cq$bfeduca==99),]


attach(jtpa.cq)

quantile(bfyrearn,c(0.33,0.67))
summary(bfyrearn)

income.braket <- function(x){
  if(x<=1000){y=1}
  else if (x>1000 & x <=4000){y=2}
  else if (x>4000){y=3}
  
  return(y)
}
jtpa.cq$preearn <-sapply(bfyrearn,income.braket,simplify = TRUE)

edu.braket <- function(x){
  if(x<=8){y=1}
  else if (x>8 & x <=10){y=2}
  else if (x>10 & x <=11){y=3}
  else if (x>11 & x <=12){y=4}
  else if (x>12 & x <=14){y=5}
  else if (x>14){y=6}
  
  return(y)
}
jtpa.cq$preedu <-sapply(bfeduca,edu.braket,simplify = TRUE)


set.seed(210010)



Y<-earnings

D  <- D 


# Controls
#X      <- "female+black+othrace+ dep+q2+q3+q4+q5+q6+agelt35+agegt54+durable+lusd+husd"         # use this for tree-based methods like forests and boosted trees
X    <- model.matrix(~(male+hsorged+black+hispanic+married+wkless13+scale(age)+scale(bfeduca)+scale(bfyrearn))^2-1,data=jtpa.cq)     # use this for rlasso etc.
W <- jtpa.cq[,((ncol(jtpa.cq)-1):ncol(jtpa.cq))]

