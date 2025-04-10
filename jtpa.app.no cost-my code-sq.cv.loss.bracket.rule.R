setwd("~...")
rm(list=ls())

jtpa.cq <- read.csv("mydata.csv",header = TRUE)[,-(1:2)]
#jtpa.cq <-jtpa.cq[-which(jtpa.cq$bfeduca==99),]


attach(jtpa.cq)

quantile(bfyrearn,c(0.33,0.67))
summary(bfyrearn)

income.braket <- function(x){
  if(x<=220){y=1}
  else if (x>220 & x <=3800){y=2}
  else if (x>3800){y=3}
  
  return(y)
}
jtpa.cq$preearn <-sapply(bfyrearn,income.braket,simplify = TRUE)

edu.braket <- function(x){
  if(x<=11){y=1}
  else if (x==12) {y=2}
  else if (x>12) {y=3}
  return(y)
}
jtpa.cq$preedu <-sapply(bfeduca,edu.braket,simplify = TRUE)


######################################################
#######Estimating the Riesz representer #########
####################################################
library(glmnet)
library("MASS")

set.seed(210010)

#number of folds for cross fitting#
K <- 5

Y<-earnings
# Treatment Indicator
D  <- D 


# Controls
#X      <- "female+black+othrace+ dep+q2+q3+q4+q5+q6+agelt35+agegt54+durable+lusd+husd"         # use this for tree-based methods like forests and boosted trees
X    <- model.matrix(~(male+hsorged+black+hispanic+married+wkless13+scale(age)+scale(bfeduca)+scale(bfyrearn))^2-1,data=jtpa.cq)     # use this for rlasso etc.
W <- jtpa.cq[,((ncol(jtpa.cq)-1):ncol(jtpa.cq))]


######################################
####calculation of debiased weight####
######################################


split     <- runif(nrow(jtpa.cq))
cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  

ATE <- as.data.frame(cbind(Y, D, W, X))
samplesize<-nrow(ATE)





policy.data <-list()


for(k in 1:K){   
  
ii  <- cvgroup == k
nii <- cvgroup != k

datatrain <- as.data.frame(ATE[nii,])
datatest <- as.data.frame(ATE[ii,])


lamda_grid <- seq(0, 5, 0.1)

J <- 10

folds <- sample(rep(1:J, length.out = nrow(datatrain)))

cv_errors_1 <- numeric(length(lamda_grid))
cv_errors_0 <- numeric(length(lamda_grid))


for (i in seq_along(lamda_grid)) {
  lamda1 <- lamda_grid[i]
  error_fold_1 <- numeric(J)
  
  lamda0 <- lamda_grid[i]
  error_fold_0 <- numeric(J)
  
  for (j in 1:J) {
    train_idx <- which(folds != j)
    val_idx <- which(folds == j)
    
    cv.train <- datatrain[train_idx,]
    cv.train1 <-    cv.train[which(cv.train$D ==1),]
    cv.train0 <-    cv.train[which(cv.train$D ==0),]
    
    cv.val <-datatrain[val_idx,]
    
    
    cvfit.r1 = cv.glmnet(x = as.matrix(cv.train1[,-(1:4)]),  y=cv.train1$Y, alpha = 1)
    
    r1.fit <- predict(cvfit.r1,newx=as.matrix(cv.train[,-(1:4)]),s = "lambda.min")
    r1.predict <- predict(cvfit.r1, newx = as.matrix(cv.val[,-(1:4)]), s = "lambda.min")
    
    cvfit.r0 = cv.glmnet(x = as.matrix(cv.train0[,-(1:4)]),  y=cv.train0$Y, alpha = 1)
    r0.fit <- predict(cvfit.r0,newx=as.matrix(cv.train[,-(1:4)]),s = "lambda.min")
    r0.predict <- predict(cvfit.r0, newx = as.matrix(cv.val[,-(1:4)]), s = "lambda.min")
    
    
    
    
    PX_train <- as.matrix(cbind(rep(1,nrow(cv.train)),cv.train[-(1:4)]))
    PX_val <- as.matrix(cbind(rep(1,nrow(cv.val)),cv.val[-(1:4)]))
    D_train <- cv.train$D
    D_val <- cv.val$D
    Y_val <-cv.val$Y
    
    
    
    m_train <-2*(r1.fit-r0.fit)
    m_val <-2*(r1.predict-r0.predict)
    
    gamma1_val <- r1.predict
    gamma0_val <- r0.predict
    
    PT_train <- PX_train * D_train
    Ghat_train <- crossprod(PT_train) / nrow(PT_train)
    Phat_train <- colMeans(apply(PX_train,2,'*',m_train))
    
    a.tilda = t(crossprod(Phat_train,ginv((Ghat_train %*% Ghat_train+lamda1*Ghat_train)) %*% Ghat_train))
    
    alpha.tilda_val = PX_val %*% a.tilda
    
    pred1 <- (alpha.tilda_val * D_val * Y_val- m_val*gamma1_val)
    
    
    PT0_train <- PX_train * (1-D_train)
    Ghat0_train <- crossprod(PT0_train) / nrow(PT0_train)
    a.tilda.0 <- t(crossprod(Phat_train,  ginv((Ghat0_train %*% Ghat0_train+lamda0*Ghat0_train)) %*% Ghat0_train)) 
    alpha.tilda.0_val = PX_val %*% a.tilda.0
    
    pred0 <- alpha.tilda.0_val * (1-D_val) * Y_val-m_val*gamma0_val
    
    
    error_fold_1[j] <- mean(pred1*pred1)
    error_fold_0[j] <-  mean(pred0*pred0)
  }
  
  cv_errors_1[i] <- mean(error_fold_1)
  cv_errors_0[i] <- mean(error_fold_0)
}


best_lamda1 <- lamda_grid[which.min(cv_errors_1)]
best_lamda0 <- lamda_grid[which.min(cv_errors_0)]



################################
################################

datatrain1 <-  as.data.frame(datatrain[which(datatrain$D ==1),])
datatrain0 <-  as.data.frame(datatrain[which(datatrain$D ==0),])


fit.r1.train= cv.glmnet(x = as.matrix(datatrain1[,-(1:4)]),  y=datatrain1$Y, alpha = 1)
r1.crossfit.train <- predict(fit.r1.train, newx = as.matrix(datatrain[,-(1:4)]), s = "lambda.min")
r1.crossfit.test <- predict(fit.r1.train, newx = as.matrix(datatest[,-(1:4)]), s = "lambda.min")

fit.r0.train= cv.glmnet(x = as.matrix(datatrain0[,-(1:4)]),  y=datatrain0$Y, alpha = 1)
r0.crossfit.train <- predict(fit.r0.train, newx = as.matrix(datatrain[,-(1:4)]), s = "lambda.min")
r0.crossfit.test <- predict(fit.r0.train, newx = as.matrix(datatest[,-(1:4)]), s = "lambda.min")


PX_cf_train <- as.matrix(cbind(rep(1,nrow(datatrain)),datatrain[-(1:4)]))
PX_cf_test <- as.matrix(cbind(rep(1,nrow(datatest)),datatest[-(1:4)]))

PT_cf_train <- apply(PX_cf_train,2,'*',datatrain$D) 

Ghat_cf_train <- as.matrix(crossprod(PT_cf_train))/nrow(PT_cf_train)

P.m1_cf <- apply(PX_cf_train,2,'*',2*(r1.crossfit.train-r0.crossfit.train))
Phat_cf <- as.vector(colMeans(P.m1_cf))


a.tilda_cf <- t(crossprod(Phat_cf,  ginv((Ghat_cf_train %*% Ghat_cf_train+best_lamda1*Ghat_cf_train)) %*% Ghat_cf_train)) #actual alpha tilda#

alpha1.predict<- PX_cf_test %*% a.tilda_cf



PT0_cf_train <- apply(PX_cf_train,2,'*',(1-datatrain$D)) 
Ghat0_cf_train <- as.matrix(crossprod(PT0_cf_train))/nrow(PT0_cf_train)

a.tilda.0_cf <- t(crossprod(Phat_cf,  ginv((Ghat0_cf_train %*% Ghat0_cf_train+best_lamda0*Ghat0_cf_train)) %*% Ghat0_cf_train)) 

alpha0.predict<- PX_cf_test %*% a.tilda.0_cf 

# Fit propensity score
cvfit.pai = cv.glmnet(
  x = as.matrix(datatrain[, -(1:4)]),
  y = datatrain$D,
  family = "binomial",
  alpha = 1
)

pai.predict <-
  predict(cvfit.pai,
          newx = as.matrix(datatest[, -(1:4)]),
          s = "lambda.min",
          type = "response")


policy.data[[k]]  <-  base::cbind(r1.crossfit.test,r0.crossfit.test, alpha1.predict,
                           alpha0.predict,pai.predict,datatest[, (1:4)])

}

ATE <- base::Reduce(base::rbind, policy.data)
colnames(ATE) <- c("fitted.r1", "fitted.r0","fitted.w1", "fitted.w0","fitted.pai", "Y", "D",
                   "preearn","preedu")

ATE$index <- as.numeric(row.names(ATE))

# Orders the rows in such that index column is in ascending order
ATE <- ATE[order(ATE$index), ]
#write.csv(ATE,file ="JTPA.processed.data.with.debiased.weight.sq.cv.loss.csv" )

#ATE <- read.csv("JTPA.processed.data.with.debiased.weight.sq.cv.loss.csv")[,-1]

ATE$preearn <- jtpa.cq$preearn 
ATE$preedu<- jtpa.cq$preedu

X.base <- model.matrix(~(male+hsorged+black+hispanic+married+wkless13+scale(age)+scale(bfeduca)+scale(bfyrearn))^2-1,data=jtpa.cq)  


attach(ATE)
####plug-in approach without debiasing or crossfitting,OLS##
reg.0 <- lm(Y ~ X.base, data = ATE, subset = (D == 0))
mu.hat.0 <- predict(reg.0, newdata = data.frame(X.base))
reg.1 <- lm(Y ~ X.base, data = ATE, subset = (D == 1))
mu.hat.1 <- predict(reg.1, newdata = data.frame(X.base))

ATE$CATEX.plug.in <-mu.hat.1-mu.hat.0


#attach(ATE)
#set.seed(120339)
#plug-in approach without debiasing or crossfitting,lasso##
#reg.0 <- cv.glmnet(x=as.matrix(X.base)[which(D==0),],y=Y[which(D==0)],alpha=1)
#mu.hat.0 <- predict(reg.0, newx= X.base, s="lambda.min")
#reg.1 <- cv.glmnet(x=as.matrix(X.base)[which(D==1),],y=Y[which(D==1)],alpha=1)
#mu.hat.1 <- predict(reg.1, newx= X.base, s="lambda.min")

#ATE$CATEX.plug.in <-mu.hat.1-mu.hat.0


##our proposed weight
ATE$weight.dr <-     as.vector(
  (ATE$`fitted.r1`-ATE$`fitted.r0`)^2
  +(ATE$D*ATE$`fitted.w1`)*(ATE$Y-ATE$`fitted.r1`)-
                                       (1-ATE$D)*ATE$`fitted.w0`*(ATE$Y-ATE$`fitted.r0`)
)


##squared weight without debiasing or cross fitting###
ATE$weight.plug.in<-     as.vector(
  (ATE$CATEX.plug.in)^2
)


##naive doubly robust weight, without using riesz representer structure###
ATE$weight.dr.naive <- as.vector(
  (ATE$`fitted.r1`-ATE$`fitted.r0`)^2
  +(ATE$`fitted.r1`-ATE$`fitted.r0`)*(2*(ATE$D/ATE$`fitted.pai`)*(ATE$Y-ATE$`fitted.r1`)-
                                        2*((1-ATE$D)/(1-ATE$`fitted.pai`))*(ATE$Y-ATE$`fitted.r0`))
)



ATE$CATEX.CF <-     as.vector(ATE$`fitted.r1` -ATE$`fitted.r0`)

ATE$CATEX.DR <-   as.vector(
  (ATE$D * ATE$Y) / ATE$`fitted.pai` - ((1-ATE$D) * ATE$Y) / (1-ATE$`fitted.pai`)-
    ( ATE$D -  ATE$`fitted.pai`)*( ATE$`fitted.r1` /ATE$`fitted.pai` + ATE$`fitted.r0`  /(1-ATE$`fitted.pai`))
  )

ATE$IPW <- as.vector(((ATE$Y*ATE$D)/(2/3))-(ATE$Y*(1-ATE$D)/(1-2/3)))



###################################
#####Decision rule################
###################################

decision.rule<- function(x,y){
  
  sq<- function(d){
    
    A1<-sum(ATE$weight.dr[which(ATE$preearn==x & ATE$preedu==y & ATE$CATEX.CF >= 0  )])
    A2<- sum(ATE$weight.dr[which(ATE$preearn==x & ATE$preedu==y)])   
    return(A2*d*d-2*A1*d) 
  }
  
  d.proposed <- round(optimise(sq,c(0,1),maximum = FALSE)$minimum,2)
  
  sq.naive<- function(d){
    
    A1<-sum(ATE$weight.dr.naive[which(ATE$preearn==x & ATE$preedu==y & ATE$CATEX.CF >= 0  )])
    A2<- sum(ATE$weight.dr.naive[which(ATE$preearn==x & ATE$preedu==y)])   
    return(A2*d*d-2*A1*d) 
  }
  
  d.dr.naive <- round(optimise(sq.naive,c(0,1),maximum = FALSE)$minimum,2)
  
  
  d.plug.in <- sum(ATE$weight.plug.in[
    which(ATE$preearn==x & ATE$preedu==y & ATE$CATEX.plug.in >= 0)])/sum(
      ATE$weight.plug.in[which(ATE$preearn==x & ATE$preedu==y)]) 
  
  cate.mean<-mean(ATE$CATEX.plug.in[which(ATE$preearn==x & ATE$preedu==y)])
  
  cate.dr.mean<-mean(ATE$CATEX.DR[which(ATE$preearn==x & ATE$preedu==y)])
  ipw.mean <-mean(ATE$IPW[which(ATE$preearn==x & ATE$preedu==y)])
  
  return(c(d.proposed,d.plug.in,cate.mean,cate.dr.mean,d.dr.naive,ipw.mean))
  
}



squared.regret.rule<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    squared.regret.rule[j,i]=decision.rule(j,i)[1]
    
  }
  
}

squared.regret.plug.in.rule<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    squared.regret.plug.in.rule[j,i]=round(decision.rule(j,i)[2],2)
    
  }
  
}


athey.wager.rule<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    athey.wager.rule[j,i]=as.numeric(decision.rule(j,i)[4]>=0)
    
  }
  
}

###using debiased structure but without using the riesz representer structure
dr.naive.rule<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    dr.naive.rule[j,i]=as.numeric(decision.rule(j,i)[5])
    
  }
  
}


CATE.DR.estimates<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    CATE.DR.estimates[j,i]=round(as.numeric(decision.rule(j,i)[4]),0)
    
  }
  
}

CATE.plug.in.estimates<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    CATE.plug.in.estimates[j,i]=round(as.numeric(decision.rule(j,i)[3]),0)
    
  }
  
}

IPW.estimates<- matrix(0,3, max(ATE$preedu)-min(ATE$preedu)+1)

for (i in 1:(max(ATE$preedu)-min(ATE$preedu)+1)) {
  for (j in 1:3) {
    IPW.estimates[j,i]=round(as.numeric(decision.rule(j,i)[6]),0)
    
  }
  
}



#####################################
##########prepare table##############
#####################################


table <- melt(squared.regret.rule)
colnames(table) <-c("pre-program income","pre-program education","squared regret rule (debiased/IV)")
table$'squared regret rule (plug-in/OLS)' <- melt(squared.regret.plug.in.rule)[,-(1:2)]
table$'CATE (plug-in)' <- round(melt(CATE.plug.in.estimates)[,-(1:2)],0)
table$'CATE (debiased)'<-round(melt(CATE.DR.estimates)[,-(1:2)],0)
table$'CATE (IPW)' <- round(melt(IPW.estimates)[,-(1:2)],0)
table

library(xtable)

table_latex <- xtable(table)
print(table_latex, type = "latex")



############################
#######heatmap plot#########
############################

library(reshape)                                                # Load reshape package

library(ggplot2)   

our.rule.plot <- melt(squared.regret.rule)
colnames(our.rule.plot) <-c("pre-program income","pre-program education","decision rule")

ggp.1 <- ggplot(our.rule.plot, aes(`pre-program education`,`pre-program income`)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = `decision rule`))
ggp.1 + scale_fill_gradient(low = "white", high = "#d12168",lim = c(0,1))+
  labs(x = "pre-program education bracket", 
       y = "pre-program income bracket", 
       title = "Our proposed approach")+
  annotate("text", x=6, y=1, label= "0.47")+
  annotate("text", x=3, y=2, label= "0.26")+
  annotate("text", x=3, y=3, label= "0.64")  


athey.wager.rule.plot <- melt(athey.wager.rule)
colnames(athey.wager.rule.plot) <-c("pre-program income","pre-program education","decision rule")

ggp.2 <- ggplot(athey.wager.rule.plot, aes(`pre-program education`,`pre-program income`)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = `decision rule`))

ggp.2 + scale_fill_gradient(low = "white", high = "#d12168",lim = c(0,1))+
  labs(x = "pre-program education bracket", 
       y = "pre-program income bracket", 
       title = "Athey-Wager approach")


squared.regret.plug.in.plot <- melt(squared.regret.plug.in.rule)
colnames(squared.regret.plug.in.plot) <-c("pre-program income","pre-program education","decision rule")

ggp.3 <- ggplot(squared.regret.plug.in.plot, aes(`pre-program education`,`pre-program income`)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = `decision rule`))

ggp.3 + scale_fill_gradient(low = "white", high = "#d12168",lim = c(0,1))+
  labs(x = "pre-program education bracket", 
       y = "pre-program income bracket", 
       title = "Squared regret plug-in approach")



dr.naive.plot <- melt(dr.naive.rule)
colnames(dr.naive.plot) <-c("pre-program income","pre-program education","decision rule")

ggp.4 <- ggplot(dr.naive.plot, aes(`pre-program education`,`pre-program income`)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = `decision rule`))

ggp.4 + scale_fill_gradient(low = "white", high = "#d12168",lim = c(0,1))+
  labs(x = "pre-program education bracket", 
       y = "pre-program income bracket", 
       title = "Debiased approach w/ naive estimation of RR")







