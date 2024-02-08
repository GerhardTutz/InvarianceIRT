


###### Programs Specific Objectivity Binary

  ###  Content:

### EstItemSepAllComp: pairwise separation
### EstItemSep: pairwise separation, only difference to first item
### EstItemPairwise: classical pairwise
### response functions
### Loglikelihood function

############################################

EstItemSepAllComp <-function(resp,I,P,gam,indsel,losspers){
  
  ##### uses all differences between pairs
  
  ### resp PxI matrix observations
  #   I number of items 
  #   P number of persons
  #   gam vector of gamma values
  #   indicator selection indsel=="q" (quadratic),indsel=="KL"(Kullback Leibler)
  #   losspers  number of persons used in minimizing loss function (first losspers persons) 
  
  respsum<-rowSums(resp)/I
  respind<- matrix(0,P,1)
  for (p in 1:P) if((respsum[p] >0) & (respsum[p] <1))respind[p,1]<-1
  
  #### Pairs omitted if  nonomit 1
  nonomit<-0
  
  if(nonomit==1){
  
  lpairs<-I*(I-1)/2
  pairs <- matrix(0,lpairs,6)
  count<- 0
  I1<-I-1
  for (i in 1:I1){i1<-i+1
  for (j in i1:I){count<-count+1
  pairs[count,1]<-i
  pairs[count,2]<-j
  }}
  
    #pairs: 1:it1 2:it2 3:n10 4:n01 5:sum 6:est
    for (c in 1:lpairs){
      n10<-0
      n01<-0
      for (r in 1:P) { if( (resp[r,pairs[c,1]]==1) & (resp[r,pairs[c,2]]==0))n10<-n10+1
      if( (resp[r,pairs[c,1]]==0) & (resp[r,pairs[c,2]]==1))n01<-n01+1   
      } 
      pairs[c,3]<-n10
      pairs[c,4]<-n01
      pairs[c,5]<-n10+n01
      pairs[c,6]<- (pairs[c,3]-pairs[c,4])/pairs[c,5]
    }
    
    pairs
  } ### end omit
  
  
  ##### Estimation: find gamma
  
  #gamtild<- c(2)
  
  lengthgam<- length(gam)
  KL<-matrix(0,lengthgam,1)
  Quadr<-matrix(0,lengthgam,1)
  #Truedis<-matrix(0,lengthgam,1)
  checkn<- matrix(0,P,I)
  
  estitemtotal<-matrix(0,lengthgam,I)
  
  ### loop gam
  
  for (l in 1:lengthgam) {
    #gammatild<-gamtild[l]
    gamma <-gam[l]
    gamma0 <-respfctinv(gamma)
    gamma1 <-respfctinv(1-gamma)
    
    
      
      ### estimation items 1 versus all
      #estitem<- c(0,pairs[1:I1,6])
      #est<- (gamma1-gamma0)*estitem
      #estitemtotal[l,]<-est
    
    
    
    
    #### all pairs
    estitemall<-matrix(0,I,I) 
    listit<-1:I
    
    for (it in 1:I){
      redit<-listit[-it]
      estitloc<-rep(0,I)
      for (c in redit) {
        sumn<-sum(abs(resp[,it]-resp[,c]))
        #sumn<-P  ### no weights!
        estitloc[c]<-(gamma1-gamma0)*(sum(resp[,it])-sum(resp[,c]))/sumn
      }
      estitemall[it,] <-  estitloc - estitloc[1] 
    } ## end it
    
    estitemtotal[l,]<-colSums(estitemall)/I
    est<-estitemtotal[l,]
    
    # person parameters
    
    P1<-P
    if (losspers > 0) P1<-losspers
    
    estpers<-matrix(0,P1,3)
    
    
    
    for (p in 1:P1){
      
      if(respind[p]==1){
      respn<-resp[p, ]
      checkn[p,]<-respn
      
      fit <- optim(0.2, Loglikper, gr = NULL, respn=respn,I=I,est=est, method = "CG",
                   lower = -Inf, upper = Inf, hessian = FALSE) 
      estpers[p,1] <-fit$par
      
      #Loglikper(4.12,respn,I,est)
      
      for (i in 1:I){
        KL[l,1] <- KL[l,1]-(resp[p,i]*log(respfctc(fit$par-est[i]))+(1-resp[p,i])*log(1-respfctc(fit$par-est[i])))
        Quadr[l,1]<-Quadr[l,1]+2*(resp[p,i]-respfctc(fit$par-est[i]))^2
      }
      
      ##checks
      
      }  ### end if respind[p]==1
      
      #for (i in 1:I){Truedis[l,1]<-Truedis[l,1]+(it[i]-est[i])^2}
      
      
    }### end p 
    #cbind(checkn,estpers[,1],pers)
    
  }### end loop find gamma
  ################################
  
  ### helpful plots 
  
  #plot(gam,KL)
  #plot(gam,Quadr)
  #plot(gam,Truedis)
  
  
  ### select gamma
  #sel<- which.min(KL)
  sel<- which.min(Quadr)
  
  if(indsel=="KL")sel<- which.min(KL)
  if(indsel=="q")sel<- which.min(Quadr)
  itfinal<-estitemtotal[sel,]
  
  #plot(it,itfinal)
  #lines(c(-10,10),c(-10,10))
  
  #cor(it,itfinal)
  
  newList <- list("itempar" = itfinal,"Quadr"=Quadr,"KL"=KL,"selectedgamma"=sel)
  return(newList)
  
  
}
################# end function




EstItemSep <-function(resp,I,P,gam,indsel){
  
  ###  uses only differences to reference item 1
  #####
  ### resp PxI matrix observations
  #   I number of items 
  #   P number of persons
  #   gam vector of gamma values
  #   indicator selection indsel=="q" (quadratic),indsel=="KL"(Kullback Leibler)
  
  
  respsum<-rowSums(resp)/I
  respind<- matrix(0,P,1)
  for (p in 1:P) if((respsum[p] >0) & (respsum[p] <1))respind[p,1]<-1
  
  #### Pairs
  
  lpairs<-I*(I-1)/2
  pairs <- matrix(0,lpairs,6)
  count<- 0
  I1<-I-1
  for (i in 1:I1){i1<-i+1
  for (j in i1:I){count<-count+1
  pairs[count,1]<-i
  pairs[count,2]<-j
  }}
  
  
  ##### Estimation: find gamma
  
  #gamtild<- c(2)
  
  lengthgam<- length(gam)
  KL<-matrix(0,lengthgam,1)
  Quadr<-matrix(0,lengthgam,1)
  Truedis<-matrix(0,lengthgam,1)
  checkn<- matrix(0,P,I)
  
  estitemtotal<-matrix(0,lengthgam,I)
  
  ### loop gam
  
  for (l in 1:lengthgam) {
    #gammatild<-gamtild[l]
    gamma <-gam[l]
    gamma0 <-respfctinv(gamma)
    gamma1 <-respfctinv(1-gamma)
    
    
    #pairs: 1:it1 2:it2 3:n10 4:n01 5:sum 6:est
    for (c in 1:lpairs){
      n10<-0
      n01<-0
      for (r in 1:P) { if( (resp[r,pairs[c,1]]==1) & (resp[r,pairs[c,2]]==0))n10<-n10+1
      if( (resp[r,pairs[c,1]]==0) & (resp[r,pairs[c,2]]==1))n01<-n01+1   
      } 
      pairs[c,3]<-n10
      pairs[c,4]<-n01
      pairs[c,5]<-n10+n01
      pairs[c,6]<- (pairs[c,3]-pairs[c,4])/pairs[c,5]
    }
    
    pairs
    
    
    ### estimation items 1 versus all
    estitem<- c(0,pairs[1:I1,6])
    est<- (gamma1-gamma0)*estitem
    estitemtotal[l,]<-est
    
    
    
    # person parameters
    
    estpers<-matrix(0,P,3)
    
    for (p in 1:P){
      
      #if(respind[p]==1){
      respn<-resp[p, ]
      checkn[p,]<-respn
      
      fit <- optim(0.2, Loglikper, gr = NULL, respn=respn,I=I,est=est, method = "CG",
                   lower = -Inf, upper = Inf, hessian = FALSE) 
      estpers[p,1] <-fit$par
      
      for (i in 1:I){
        KL[l,1] <- KL[l,1]-(resp[p,i]*log(respfctc(fit$par-est[i]))+(1-resp[p,i])*log(1-respfctc(fit$par-est[i])))
        Quadr[l,1]<-Quadr[l,1]+2*(resp[p,i]-respfctc(fit$par-est[i]))^2
      }
      ##checks
      
      #}  ### end if respind[p]==1
      
      #for (i in 1:I){Truedis[l,1]<-Truedis[l,1]+(it[i]-est[i])^2}
      
      
    }### end p 
    #cbind(checkn,estpers[,1],pers)
    
  }### end loop find gamma
  ################################
  
  ### helpful plots 
  
  #plot(gam,KL)
  #plot(gam,Quadr)
  #plot(gam,Truedis)
  
  
  ### select gamma
  #sel<- which.min(KL)
  sel<- which.min(Quadr)
  
  if(indsel=="KL")sel<- which.min(KL)
  if(indsel=="q")sel<- which.min(Quadr)
  itfinal<-estitemtotal[sel,]
  
  #plot(it,itfinal)
  #lines(c(-10,10),c(-10,10))
  
  #cor(it,itfinal)
  
  newList <- list("itempar" = itfinal,"Quadr"=Quadr,"KL"=KL,"selectedgamma"=sel)
  return(newList)
  
  
}
################# end function


EstItemPairwise <-function(resp,I,P){
  
  ##### uses all differences between pairs
  
  ### resp PxI matrix observations
  #   I number of items 
  #   P number of persons
  #   gam vector of gamma values
  #   indicator selection indsel=="q" (quadratic),indsel=="KL"(Kullback Leibler)
  
  
  gam<-c(.1)  #  fixed, since not needed
  
  respsum<-rowSums(resp)/I
  respind<- matrix(0,P,1)
  for (p in 1:P) if((respsum[p] >0) & (respsum[p] <1))respind[p,1]<-1
  
  #### Pairs
  
  lpairs<-I*(I-1)/2
  pairs <- matrix(0,lpairs,6)
  count<- 0
  I1<-I-1
  for (i in 1:I1){i1<-i+1
  for (j in i1:I){count<-count+1
  pairs[count,1]<-i
  pairs[count,2]<-j
  }}
  
  
  ##### Estimation: find gamma
  
  #gamtild<- c(2)
  
  lengthgam<- length(gam)
  KL<-matrix(0,lengthgam,1)
  Quadr<-matrix(0,lengthgam,1)
  #Truedis<-matrix(0,lengthgam,1)
  checkn<- matrix(0,P,I)
  
  estitemtotal<-matrix(0,lengthgam,I)
  
  ### loop gam
  
  for (l in 1:lengthgam) {
    #gammatild<-gamtild[l]
    gamma <-gam[l]
    gamma0 <-respfctinv(gamma)
    gamma1 <-respfctinv(1-gamma)
    
    
    
    #### alternative
    estitemall<-matrix(0,I,I) 
    listit<-1:I
    
    
    #### now just once 
    I<-1
    for (it in 1:I){
      redit<-listit[-it]
      estitloc<-rep(0,I)
      for (c in redit) {
        sumn<-sum(abs(resp[,it]-resp[,c]))
        inddiff<-abs(resp[,it]-resp[,c])
        
        #estitloc[c]<-(gamma1-gamma0)*(sum(resp[,it])-sum(resp[,c]))/sumn
        ### here
        #estitloc[c]<-log(resp[,it]%*%inddiff)-log(resp[,c]%*%inddiff)
        estitloc[c]<-log((resp[,it]%*%inddiff+.01)/(resp[,c]%*%inddiff+.01))
      }
      estitemall[it,] <-  estitloc - estitloc[1] 
    } ## end it
    
    estitemtotal[l,]<-colSums(estitemall)/I
    
    # person parameters
    
    
    
    
  }### end loop find gamma
  ################################
  
  
  itfinal<-estitemtotal[1,]
  
  newList <- list("itempar" = itfinal)
  return(newList)
  
  
}
################# end function







###### response fctions
respfctcNV <-function(x){
  r<-pnorm(x) 
  return(r)  
}
respfctinvNV <-function(x){
  r<-qnorm(x) 
  return(r)
}

respfctclogit <-function(x){
  r<-exp(x)/(1+exp(x)) 
  return(r)  
}
respfctinvlogit <-function(x){
  r<-log(x/(1-x)) 
  return(r)  
}

respfctcGumbel <-function(x){
  r<-exp(-exp(-x)) 
  return(r)  
}
respfctinvGumbel <-function(x){
  r<--log(-log(x)) 
  return(r)  
}


respfctcGompertz <-function(x){
  r<-1-exp(-exp(x)) 
  return(r)  
}

respfctinvGompertz <-function(x){
  r<-log(-log(1-x)) 
  return(r)  
}


#### lloglik function

Loglikper <-function(par,respn,I,est){
  
  ridge<-0.05
  c0<-0.001  ## correct for existence
  sum<-0
  for(i in 1:I){
    respfctcorr<-c0+(1-2*c0)*respfctc(par-est[i])
    sum<- sum+ respn[i]*log(respfctcorr)+(1-respn[i])*log(1-respfctcorr)  
    sum<- sum - ridge*par^2  ##ridge
  }
  sum<--sum
  
  return(sum)}
#### end function




