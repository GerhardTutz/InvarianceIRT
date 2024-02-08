

#### Simulation Spec Obj with function

library("eRm")
library(spatstat)  ### gauss Hermite

## run first:

source("Responsefunctions.R")
source("ProgramsSepEstimation.R")
source("ProgramsPoly2.R")


### items and persons
I <- 6   # number of items

it<-c(0, -1.5, -1,   0.5, 1.2, 1.5)  # item parameters
itm<-matrix(it)
K<-1

P<-200   # number of persons

#### alternative items 1

it<-c(it,it,it)
I<-length(it)
P<-100


#### alternative items 2

Igen<- 17
itn<-0
for (i in 1:Igen)
itn<- c(itn,2*sin(-1.5+.25*(i-1)))
plot(seq(1,Igen+1,1),itn)

it<-itn
I<-Igen+1

P<-800
######################



#### choose

gam<- seq(.03,.4,.015 )  ### gamma values

indsel<-"KL"  ## choice selection


##### choose response function

respf<- "NV"
#respf<- "logistic"
#respf<- "Gumbel"
#respf<- "Gompertz"



if(respf=="NV"){respfctd <- respfctdNV
respfctc<-respfctcNV
respfctder<-respfctderNV
respfctinv<-respfctinvNV}

if(respf=="logistic"){respfctd <- respfctdlogit
respfctc<-respfctclogit
respfctder<-respfctderlogit
respfctinv<-respfctinvlogit}

if(respf=="Gumbel"){respfctd <- respfctdGumbel
respfctc<-respfctcGumbel
respfctder<-respfctderGumbel
respfctinv<-respfctinvGumbel}

if(respf=="Gompertz"){respfctd <- respfctdGompertz
respfctc<-respfctcGompertz
respfctder<-respfctderGompertz
respfctinv<-respfctinvGompertz}






####simulation

numsim<-5

RaschInd<-1  ### if 1 fits conditional Rasch
marginalInd<-0   ### if 1 fits marginal model
poly<-1  ### if 1 fits several times with changing anchor items


itsim<-matrix(0,numsim,I)
unconditsim<-matrix(0,numsim,I)   ### conditional estimate
itsimpair<-matrix(0,numsim,I)
itsimmarginal<-matrix(0,numsim,I)


##1: separate all 2: all with poly 3:conditional Rasch 4:pairwise 5:marginal
Abserr<-matrix(0,numsim,5)
 
#AbserrPart<-matrix(0,numsim,4)  ### compares versions of part

failcond<-0
numsimcond<-0

for (s in 1:numsim){
  
  ##### data generation
  pers <- rnorm(P, mean = 0, sd = 1)   #### person par
  
  #pers<-pers^2  #### squared
  resp<- matrix(0,P,I)
  
  for (i in 1:I){
    for (p in 1:P){
      prob <- respfctc(pers[p]-it[i]) 
      draw<- runif(1, min = 0, max = 1)
      if(draw < prob)resp[p,i]<-1
    }}
  
  
  ### indicator for total scores
  respsum<-rowSums(resp)/I
  respind<- matrix(0,P,1)
  for (p in 1:P) if((respsum[p] >0) & (respsum[p] <1))respind[p,1]<-1
  
  
  ## separation estimator
  Trial<-EstItemSepAllComp(resp,I,P,gam,indsel,losspers=P)
  itsim[s,] <- Trial$itempar
  Abserr[s,1]<-sum(abs(Trial$itempar-it))/I
  
  
  
  if(poly==1){TrialPoly<-EstItemSepPolyMultiLossAll(resp,I,P,K,gam,indsel,weight="w",multismooth=0)
  Abserr[s,2] <-   sum(abs(TrialPoly$itempar-it))/I
  itfinal<-TrialPoly$itempar
  itsim[s,]<-itfinal}
  
  
  if(marginalInd==1) {Trialmarg<-MargFit(resp,I)
  itsimmarginal[s,]<- Trialmarg$itempar
  Abserr[s,5] <-   sum(abs(Trialmarg$itempar-it))/I
  }
    
   
  ############################
  
  
  #### plot for true values - takes time, only for check=1
  check<-0
  if (check==1){
  Truedis<-matrix(0,length(gam),1)
  
  for (l in 1:length(gam)) {
    Trialn<-EstItemSep(resp,I,P,gam[l],indsel=P)
    for (i in 1:I){Truedis[l,1]<-Truedis[l,1]+(it[i]-Trialn$itempar[i])^2}
               }
  plot(gam,Truedis)
    }
  
  
  #### pairwise conditional
  
  if (RaschInd==1) {Trialpair<-EstItemPairwise(resp,I,P)
     itsimpair[s,]<-Trialpair$itempar}

  
  #### conditional 
  
  if(RaschInd ==1){res.rasch <- RM(resp,sum0="FALSE")
    I1<-I-1
  if(length(res.rasch$etapar) == I1) {numsimcond<-numsimcond +1
    unconditsim[numsimcond,]<-c(0,res.rasch$etapar)
    Abserr[numsimcond,3] <-   sum(abs(unconditsim[s,]-it))/I 
                                     }
    
    if(length(res.rasch$etapar) < I1) failcond<-failcond+1
    }
  
  #pres.rasch <- person.parameter(res.rasch)
  
  Abserr[s,4] <-   sum(abs(itsimpair[s,]-it))/I
  
  print(s)
  
  } ### end loop sim

############################
####################################### end simulation






### summary simulation



##1: separate all 2: all with poly 3:conditional Rasch 4:pairwise 5:marginal

Error<-colSums(Abserr)/numsim
Error[3]<-sum(Abserr[1:numsimcond,3])/numsimcond  ### only valid for conditional
Error

mem<-colSums(Abserr)/numsim



failcond  ### number of simulations where conditional failed!!!
numsimcond  ### number valid estimates conditional

########################################################




## std err 

stderr<-matrix(0,1,I)
for (i in 1:I){stderr[1,i]<-sd(itsim[,i])}
stderr

#############
###########################################





###   plots loss functions

plot(gam,Trial$KL,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="Kullback-Leibler",xlab ="gamma", type="b",
     ylab="")



### compare loss functions with different numbers of persons (if estimated)

plot(gam,Trial$KL/P,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="Kullback-Leibler",xlab ="gamma", type="b",
     ylab="")
#lines(gam,Trialredpers100$KL/100, type="b")
#lines(gam,Trialredpers50$KL/50, type="b")

plot(gam,Trial$Quadr/P,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="Quadratic",xlab ="gamma", type="b",
     ylab="")
#lines(gam,Trialredpers100$Quadr/100, type="b")
#lines(gam,Trialredpers50$Quadr/50, type="b")

Trial$selectedgamma
#Trialredpers100$selectedgamma
#Trialredpers50$selectedgamma

#### check true loss

Truedis<-matrix(0,length(gam),1)
for (l in 1:length(gam)) {
  Trialn<-EstItemSep(resp,I,P,gam[l],indsel=P)
  for (i in 1:I){Truedis[l,1]<-Truedis[l,1]+(it[i]-Trialn$itempar[i])^2}
}
plot(gam,Truedis)




##### boxplots simulations


## separation
dum<-seq(1,I,1)
boxplot(itsim,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Pairwise separation estimates",xlab ="Items", type="b",
        ylab="") #,names=c("1", "2", "3","4","5","6"))
lines(dum,it,pch=16,type='p',cex=1.3)

## conditional
boxplot(unconditsim[1:numsimcond,],cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Conditional estimates",xlab ="Items", type="b",
        ylab="")#,names=c("1", "2", "3","4","5","6"))
lines(dum,it,pch=16,type='p',cex=1.3)

# pairwise
boxplot(itsimpair,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Pairwise conditional estimates",xlab ="Items", type="b",
        ylab="")#,names=c("1", "2", "3","4","5","6"))
lines(dum,it,pch=16,type='p',cex=1.3)

# marginal
boxplot(itsimmarginal,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Marginal estimates",xlab ="Items", type="b",
        ylab="")#,names=c("1", "2", "3","4","5","6"))
lines(dum,it,pch=16,type='p',cex=1.3)


## without main
dum<-seq(1,I,1)
boxplot(itsim,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5,xlab ="Items", type="b",
        ylab="") #,names=c("1", "2", "3","4","5","6"))
lines(dum,it,pch=16,type='p',cex=1.3)



### with ylim fixed

ylims<-c(-2.2,2.8)
boxplot(itsim,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Pairwise separation estimates",xlab ="Items", type="b",
        ylab="", ylim=ylims)
lines(dum,it,pch=16,type='p',cex=1.3)

boxplot(unconditsim,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Conditional estimates",xlab ="Items", type="b",
        ylab="", ylim=ylims)
lines(dum,it,pch=16,type='p',cex=1.3, ylim=ylims)

# pairwise
boxplot(itsimpair,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="Pairwise conditional estimates",xlab ="Items", type="b",
        ylab="", ylim=ylims)
lines(dum,it,pch=16,type='p',cex=1.3, ylim=ylims)



##### plots density items
itplot<-5
plotit<-itsim[,itplot]

d <- density(plotit, bw = "sj",adjust = 1.4)
xlims<-c(0,2.8)
plot(d,main ="",cex=1.3,xlab ="",cex.lab=1.5,cex.axis=1.5,xlim=xlims)
lines(it[itplot],.02,pch=16,type="p")









#################################################    
##  Bootstrap
numboot<-50
estboot<-matrix(0,numboot,I)
rep<-200
stbootrep<- matrix(0,rep,I)

for (r in 1:rep){
 for (b in 1:numboot){
 indpers<-seq(1,P,1)
 resample <- sample(indpers, replace = TRUE)
 resample

 newresp <- resp[resample,] 

 Trialb<-EstItemSepAllComp(newresp,I,P,gam,indsel="KL",losspers=P)
 estboot[b,]<-Trialb$itempar
}
summary(estboot)
stdboot<-matrix(0,1,I)
for (i in 1:I){stdboot[1,i]<-sd(estboot[,i])}
stbootrep[r,]<-stdboot
} ## end rep



stdboot  
stderr  


#### plots

plot(seq(2,I,1),stderr[2:I],pch=16,type="p")
for(r in 1:rep)lines(seq(2,I,1),stbootrep[r,2:I],type="p")


boxplot(stbootrep[,2],stbootrep[,3],stbootrep[,4],stbootrep[,5],stbootrep[,6],cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="",xlab ="Items", type="b",
        ylab="",names=c( "2", "3","4","5","6"))
lines(seq(2,I,1),stderr[2:I],pch=16,type="p",cex=1.3)


stbootrepred<-stbootrep[,2:I]
stderrred<-stderr[2:I]
boxplot(stbootrepred,cex.axis=1.5,cex.lab=1.5,
        cex.main=1.5, main ="",xlab ="Items", type="b",
        ylab="",names=c( "2", "3","4","5","6"),ylim=c(.12,.55))
lines(seq(1,I-1,1),stderrred,pch=16,type="p",cex=1.3)
colMeans(stbootrepred)






### plots approximate model


xpl<-seq(-3,3,0.1)
ypl<-seq(-1,1,0.1)
for (r in 1:length(xpl))ypl[r]<-respfctcGumbel(xpl[r]) 
plot(xpl,ypl)

gamma<-0.05 
theta<-seq(-2,2,.1)
thetastar<-0*theta
for (i in 1:length(theta)) thetastar[i]<-respfctinv(respfctc(theta[i])*(1-gamma)+gamma) 

plot(theta,thetastar,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="",xlab ="theta",ylab ="theta*", type="b")
lines(c(-3,3),c(-3,3))


##delta
gamma<-0.05 #fixed
thetan<-1  #fixed
thetastarn<-respfctinv(respfctc(thetan)*(1-gamma)+gamma)

delt<-seq(-1.5,1.5,0.1)
deltstar<- 0*delt 


for (i in 1:length(delt)) deltstar[i]<--respfctinv(respfctc(thetan-delt[i])*(1-gamma)+gamma)+ thetastarn

plot(delt,deltstar,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="",xlab ="delta",ylab ="delta*", type="b")
lines(c(-3,3),c(-3,3))


plot(delt,deltstar-delt,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="",xlab ="delta",ylab ="delta*", type="b")
lines(c(-3,3),c(0,0))

#### old
Trial$KL-Trialredpers$KL


plot(gam,Trial$Quadr,cex.axis=1.5,cex.lab=1.5,
     cex.main=1.5, main ="Quadratic",xlab ="gamma", type="b",
     ylab="")


