

## Programs polynomous responses (cum and seq) Rest possibly in ProgramsPoly.R



###################################

EstItemSepPolyMultiLoss <-function(resp,I,P,K,gam,indsel,weight,multismooth){
  
  ### cumulative model with multivariate loss  
  ###  uses only differences to reference item 1
  #####
  ### resp PxI matrix observations
  #   I number of items 
  #   P number of persons
  #   gam vector of gamma values
  #   indicator selection indsel=="q" (quadratic),indsel=="KL"(Kullback Leibler)
  #   if weight=="w" variable weights (denominator) )
  #   if multismooth > 0 smoothing parameter adapted (2*multismooth neighbourhood values)
  
  
  respsum<-rowSums(resp)
  respind<- matrix(0,P,1)
  for (p in 1:P) if((respsum[p] >0) & (respsum[p] <I*K))respind[p,1]<-1
  
  lengthgam<-length(gam) 
  
  #### Pairs responsecat
  #
  lpairs<-(K*I)
  #lpairs<-(K*(K-1)/2)*((I*(I-1)/2)+I)
  pairs <- matrix(0,lpairs,6)  ## in 5 Diff
  
  ## for gam
  
  lengthgam<- length(gam)
  KL<-matrix(0,lengthgam,1)
  Quadr<-matrix(0,lengthgam,1)
  Truedis<-matrix(0,lengthgam,1)
  checkn<- matrix(0,P,I)
  
  estitemtotal<-array(0,dim=c(I,K,lengthgam))
  
  #### pairs
  Kr<-K-1
  Ir<-I-1
  countpairs<-0
  
  
  ##### in the Following i=r=1 !
  for (r in 1:1){
    r1<-r+1
    for (q in 1:K){
      
      for (j in 1:I){
        
        #for (i in 1:I){ here
        for (i in 1:1){ 
          countpairs<-countpairs+1
          pairs[countpairs,1]<-r
          pairs[countpairs,2]<-q
          pairs[countpairs,3]<-j
          pairs[countpairs,4]<-i  
          
          count<-0
          count1<-0
          
          for (p in 1:P)  {if((resp[p,i] >= r)&(resp[p,j] < q)) count<-count+1
          if((resp[p,i] < r)&(resp[p,j] >= q))count<-count-1
          if(((resp[p,i] >= r)&(resp[p,j] < q)) |(resp[p,i] < r)&(resp[p,j] >= q))count1<-count1+1} ###!!
          
          pairs[countpairs,5]<-count
          pairs[countpairs,6]<-count1
          #if((i==j)&(r >=q))pairs[countp,5]<-0  ## impossible case
        }}
    }}
  
  pairs
  norm<-max(abs(pairs[,5]))  ### normierung
  
  
  
  ###############################  
  ## gam selection
  estpers<-matrix(0,P,1)
  
  for (l in 1:lengthgam) {
    
    gamma <-gam[l]
    gamma0 <-respfctinv(gamma)
    gamma1 <-respfctinv(1-gamma)
    
    
    ### estimation items 1 versus all
    
    for (j in 1:I){
      for (q in 1:K){
        
        countnum<-0
        
        est<-0
        for (s in 1:countpairs){if((pairs[s,1]==1)&(pairs[s,2]==q)&(pairs[s,3]==j)&(pairs[s,4]==1))
          est<-est+pairs[s,5] 
        countnum<-countnum+abs(pairs[s,5])
        norm1<-pairs[s,6]}
        
        
        if(countnum >0)estitemtotal[j,q,l]<-(gamma1-gamma0)*est/norm  #countnum    
        if((countnum >0) & (weight=="w"))estitemtotal[j,q,l]<-(gamma1-gamma0)*est/norm1 
      }}
    #max(abs(pairs[,5]))
    #estitemtotal[1,1,l]<-0 
    arg<-1:K
    #  monotony
    for(ii in 1:I) {rr<-isoreg(arg,estitemtotal[ii,,l])
    estitemtotal[ii,,l]<- rr$yf+arg*0.001}
    
    #estitemtotal[,,l]
    
    
    for (p in 1:P){
      #if(respind[p]==1){
      
      respn<-resp[p, ]
      estitnow<-estitemtotal[,,l]
      
      if(K==1){estitnow<-as.matrix(estitnow)} 
      fit <- optim(0.2, LoglikperPoly, gr = NULL, respn=respn,I=I,K=K,estitnow=estitnow, method = "CG",
                   lower = -Inf, upper = Inf, hessian = FALSE) 
      
      estpers[p,1] <-fit$par
      
      
      for (i in 1:I){
        prob0<-ProbCum(fit$par, estitemtotal[i,,l],K)
        KL[l,1] <- KL[l,1]-log(prob0[resp[p,i]+1]+.0001)
        
        respn<-resp[p,i]+1 ## in 1...K+1
        sumn<-(1-prob0[respn])^2
        listp<-c(1:K,K+1)
        redit<-listp[-respn]
        
        for(d in redit)sumn<-sumn+(prob0[d])^2
        Quadr[l,1]<-Quadr[l,1]+sumn
      }
      ##checks
      
      #}  ### end if respind[p]==1
      
    }### end person parameters 
    
  }###  gam selection
  
  
  ################################
  
  ### helpful plots 
  
  #plot(gam,KL)
  #plot(gam,Quadr)
  #plot(gam,Truedis)
  
  
  
  ### select gamma
  sel<- 1  ### dummy
  sel<- which.min(Quadr)
  crit<- Quadr[sel,1]
  if(indsel=="KL"){sel<- which.min(KL)
  crit<- KL[sel,1]}  
  
  itfinal<-estitemtotal[,,sel]
  
  
  ### multidimensional 
  
  if(multismooth > 0){ 
    
    #choice<-c(sel)
    #if ((sel > 2) & (sel <length(gam)-2))choice<-c(sel-2,sel-1,sel+1,sel+2)
    
    neighb<-c((sel-multismooth):(sel-1),(sel+1):(sel+multismooth))
    choice<-neighb
    if(min(neighb)< 1)choice<-neighb+min(neighb)+1
    if(max(neighb)> length(gam))choice<-neighb-(max(neighb)-length(gam))
    
    for (it in 1:I){
      for (c in choice){
        
        KLn<-0
        Quadrn<-0
        
        estitnow<-itfinal
        estitnow[it,]<-estitemtotal[it,,c]
        
        if(K==1)estitnow<-as.matrix(estitnow)
        
        
        for (p in 1:P){
          respn<-resp[p, ]
          fit <- optim(0.2, LoglikperPoly, gr = NULL, respn=respn,I=I,K=K,estitnow=estitnow, method = "CG",
                       lower = -Inf, upper = Inf, hessian = FALSE) 
          estpers[p,1] <-fit$par
          
          
          for (i in 1:I){
            prob0<-ProbCum(fit$par, estitnow[i,],K)
            KLn <- KLn-log(prob0[resp[p,i]+1]+.0001)
            
            respn<-resp[p,i]+1 ## in 1...K+1
            sumn<-(1-prob0[respn])^2
            listp<-c(1:K,K+1)
            redit<-listp[-respn]
            
            for(d in redit)sumn<-sumn+(prob0[d])^2
            Quadrn<-Quadrn+sumn
          } } #end p
        
        if((indsel=="KL")&(KLn <crit))itfinal<-estitnow
        if((indsel=="q")&(Quadrn <crit))itfinal<-estitnow
        
        
      } } ### end it
    
  } ### end if multismooth
  
  
  
  if (multismooth > 1)newList <- list("itempar" = itfinal,"Quadr"=Quadr,"KL"=KL,"selectedgamma"=sel, 
                                      "KLFinal (if KL selected)"=KLn,"quadrFinal (if Quadr selected)"=Quadrn)
  
  
  
  if (multismooth <= 0)newList <- list("itempar" = itfinal,"Quadr"=Quadr,"KL"=KL,"selectedgamma"=sel, 
                                       "KLFinal"=KL[sel,1])
  
  return(newList)
  
  
}
################# end function


########################
EstItemSepPolyMultiLossAll<-function(resp,I,P,K,gam,indsel,weight,multismooth){
  
  #### averaged over all anchor items
  
  estall<- matrix(0,I,K)
  
  for (i in 2:I) {respnew<-resp
  respnew[,i]<-resp[,1]
  respnew[,1]<-resp[,i]
  fit<-EstItemSepPolyMultiLoss(respnew,I,P,K,gam,indsel=indsel,weight=weight,multismooth=multismooth)
  estnew<-fit$itempar
  if(K >1){estnew[i,]<-fit$itempar[1,]
  estnew[1,]<-fit$itempar[i,]
  estnew<- estnew- estnew[1,1]}
  if(K ==1){estnew[i]<-fit$itempar[1]
  estnew[1]<-fit$itempar[i]
  estnew<- estnew- estnew[1]}
  
  
  estall<- estall+estnew/(I-1)
  }## end i
  
  newList<-list("itempar"=estall)
  
  return(newList)}
######################






################################

LoglikperPoly <-function(par,respn,I,K,estitnow){
  
  ridge<-0.05
  c0<-0.001  ## correct for existence
  sum<-0
  for(i in 1:I){
    
    probgr<-rep(0,K)
    prob<-rep(0,K)
    
    
    for (r in 1:K){probgr[r]<-respfctc(par-estitnow[i,r])}
    Kdum<-K-1
    for (r in 1:Kdum){prob[r]<-probgr[r]-probgr[r+1]}
    prob[K]<-probgr[K]
    prob0<-c(1-probgr[1],prob) #with resp o at first place
    
    #cat<-which.max(rmultinom(1, 1, prob0))
    
    #draw<- runif(1, min = 0, max = 1)
    
    
    
    sum<- sum+ log(prob0[respn[i]+1]+0.0001) 
    sum<- sum - ridge*par^2  ##ridge
  }
  sum<--sum
  
  return(sum)}
#### end function




###########################

ProbCum <-function(theta, delta,K){
  
  ### computes probabilities of cumulative model 
  ###  k+1 categories 0,1,  K
  ### delta: k-vector of thresholds
  
  probgr<-rep(0,K)
  prob<-rep(0,K)
  
  
  for (r in 1:K){probgr[r]<-respfctc(theta-delta[r])}
  Kdum<-K-1
  for (r in 1:Kdum){prob[r]<-probgr[r]-probgr[r+1]}
  prob[K]<-probgr[K]
  prob0<-c(1-probgr[1],prob) #with resp o at first place
  
  
  return(prob0)}
#### end function





###############################################
EstItemSepSeq2 <-function(resp,I,P,K,gam,indsel,weight,multismooth){
  
  ### sequential model 2 with multismooth
  ###  uses only differences to reference item 1
  #####
  ### resp PxI matrix observations
  #   I number of items 
  #   P number of persons
  #   gam vector of gamma values
  #   indicator selection indsel=="q" (quadratic),indsel=="KL"(Kullback Leibler)
  
  
  respsum<-rowSums(resp)
  respind<- matrix(0,P,1)
  for (p in 1:P) if((respsum[p] >0) & (respsum[p] <I*K))respind[p,1]<-1
  
  lengthgam<-length(gam) 
  
  #### Pairs responsecat for items larger than 1
  
  # in 3 differences in 5 weighted differences
  
  ## for gam
  
  lengthgam<- length(gam)
  KL<-matrix(0,lengthgam,1)
  Quadr<-matrix(0,lengthgam,1)
  Truedis<-matrix(0,lengthgam,1)
  checkn<- matrix(0,P,I)
  
  estitemtotal<-array(0,dim=c(I,K,lengthgam))
  
  #### pairs
  estprel<-matrix(0,I,K)
  
  lpairs<-(I-1)*K
  pairs <- matrix(0,lpairs,5)  ## in 3 Diff
  
  countpairs<-0
  
  for (q in 1:K){
    
    for (j in 2:I){
      countpairs<-countpairs+1
      pairs[countpairs,1]<-j
      pairs[countpairs,2]<-q
      
      count<-0
      countabs<-0
      for (p in 1:P)  {if((resp[p,1] >= 1)&(resp[p,j] >= q-1)) count<-count+1
      if((resp[p,j] >= q))count<-count-1
      if((resp[p,j] >= q-1))countabs<-countabs+1}
      pairs[countpairs,3]<-count
      pairs[countpairs,4]<-countabs
      pairs[countpairs,5]<-count*countabs
      estprel[j,q]<-count
      if(weight=="w")estprel[j,q]<-count/countabs #weighted
    }
  }
  
  pairs
  
  itm
  estprel
  
  #pairs2 new for item 1 with item 2
  
  lpairs2<-K-1
  pairs2 <- matrix(0,lpairs2,5)  ## in 3 Diff
  
  s<-2
  countpairs2<-0
  for (q in 2:K){
    countpairs2<-countpairs2+1
    pairs2[countpairs2,1]<-1
    pairs2[countpairs2,2]<-q
    
    count<-0
    countabs<-0
    for (p in 1:P)  {if((resp[p,2] >= 1)&(resp[p,1] >= q-1)) count<-count+1
    if(resp[p,1] >= q) count<-count-1
    
    if(resp[p,1] >= q-1)countabs<-countabs+1}
    pairs2[countpairs2,3]<-count
    pairs2[countpairs2,4]<-countabs
    pairs2[countpairs2,5]<-count/countabs
    
    estprel[1,q]<- count + estprel[2,1]  ### reference item 2
    if(weight=="w")estprel[1,q]<-count/countabs +estprel[2,1] ## weighted
  }
  
  pairs2 
  itm
  estprel
  
  
  maxpairs<-max(abs(pairs[,3]))  ### normierung  
  if(weight=="w")maxpairs<-1
  
  ###############################  
  ## gam selection
  
  estpers<-matrix(0,P,3)
  
  for (l in 1:lengthgam) {
    
    gamma <-gam[l]
    gamma0 <-respfctinv(gamma)
    gamma1 <-respfctinv(1-gamma)
    
    
    ### estimation items 1 versus all
    countdum<-0
    estitemtotal[,,l]<-(gamma1-gamma0)*estprel/maxpairs
    
    for (p in 1:P){
      #if(respind[p]==1){
      
      respn<-resp[p, ]
      estitnow<-estitemtotal[,,l]
      
      if(K==1){estitnow<-as.matrix(estitnow)} 
      fit <- optim(0.2, LoglikperSeq, gr = NULL, respn=respn,I=I,K=K,estitnow=estitnow, method = "CG",
                   lower = -Inf, upper = Inf, hessian = FALSE) 
      
      #estpers[p,l] <-fit$par  
      #LoglikperSeq(.2,resp[p, ],I,K,estitnow)
      
      #persd<-seq(-1,1,.1)
      #Logd<-matrix(0,length(persd),1)
      #for (id in 1:length(persd)) Logd[id,1]<-LoglikperSeq(persd[id],resp[p, ],I,K,estitnow)
      #plot(persd,Logd)
      
      for (i in 1:I){
        prob0<-ProbSeq(fit$par, estitemtotal[i,,l],K)
        KL[l,1] <- KL[l,1]-log(prob0[resp[p,i]+1]+.0001)
        
        respn<-resp[p,i]+1 ## in 1...K+1
        sumn<-(1-prob0[respn])^2
        listp<-c(1:K,K+1)
        redit<-listp[-respn]
        
        for(d in redit)sumn<-sumn+(prob0[d])^2
        Quadr[l,1]<-Quadr[l,1]+sumn
      }
      
    }  ### end p person parameters 
    
    
  }  ###  gam selection
  
  
  ################################
  
  ### helpful plots 
  
  #plot(gam,KL)
  #plot(gam,Quadr)
  #plot(gam,Truedis)
  
  
  ### select gamma
  sel<- 1  ### dummy
  sel<- which.min(Quadr)
  crit<- Quadr[sel,1]
  if(indsel=="KL"){sel<- which.min(KL)
  crit<- KL[sel,1]} 
  
  itfinal<-estitemtotal[,,sel]
  
  
  ### multidimensional 
  
  if(multismooth > 0){ 
    
    neighb<-c((sel-multismooth):(sel-1),(sel+1):(sel+multismooth))
    choice<-neighb
    if(min(neighb)< 1)choice<- 1:(2*multismooth)     
    if(max(neighb)> length(gam))choice<-neighb-(max(neighb)-length(gam))
    
    for (it in 1:I){
      for (c in choice){
        
        KLn<-0
        Quadrn<-0
        
        estitnow<-itfinal
        estitnow[it,]<-estitemtotal[it,,c]
        
        if(K==1)estitnow<-as.matrix(estitnow)
        
        for (p in 1:P){
          respn<-resp[p, ]
          fit <- optim(0.2, LoglikperSeq, gr = NULL, respn=respn,I=I,K=K,estitnow=estitnow, method = "CG",
                       lower = -Inf, upper = Inf, hessian = FALSE) 
          estpers[p,1] <-fit$par
          
          
          for (i in 1:I){
            prob0<-ProbSeq(fit$par, estitnow[i,],K)
            KLn <- KLn-log(prob0[resp[p,i]+1]+.0001)
            
            respn<-resp[p,i]+1 ## in 1...K+1
            sumn<-(1-prob0[respn])^2
            listp<-c(1:K,K+1)
            redit<-listp[-respn]
            
            for(d in redit)sumn<-sumn+(prob0[d])^2
            Quadrn<-Quadrn+sumn
          } } #end p
        
        if((indsel=="KL")&(KLn <crit))itfinal<-estitnow
        if((indsel=="q")&(Quadrn <crit))itfinal<-estitnow
        
        
      } } ### end it
    
  } ### end if multismooth
  
  
  
  if (multismooth == 0) newList <- list("itempar" = itfinal,"Quadr"=Quadr,"KL"=KL,"selectedgamma"=sel)
  
  
  if (multismooth > 0) newList <- list("itempar" = itfinal,"Quadr"=Quadr,"KL"=KL,"selectedgamma"=sel, 
                                       "KLFinal (if KL selected)"=KLn,"quadrFinal (if Quadr selected)"=Quadrn)
  return(newList)
  
}
################# end function
#cutoff



LoglikperSeq <-function(par,respn,I,K,estitnow){
  
  ridge<-0.05
  c0<-0.001  ## correct for existence
  
  sum<-0
  for(i in 1:I){
    
    probF<-rep(0,K)
    prob0<-rep(0,K+1)
    
    ### seq
    for (r in 1:K)probF[r]<-respfctc(par-estitnow[i,r])
    
    ##comp prob
    prob0[1]<-1-probF[1]
    
    Kn1<-K-1
    for (r in 1:Kn1){dum<-1-probF[r+1]  
    for (j in 1:r) {dum<-dum*probF[j]}
    prob0[r+1]<-dum} #with resp o at first place
    
    dumn<-1
    for (j in 1:K) dumn<-dumn*probF[j]
    prob0[K+1]<-dumn
    #sum(prob0)
    
    
    ######
    
    
    sum<- sum+ log(prob0[respn[i]+1]+0.0001) 
    sum<- sum - ridge*par^2  ##ridge
  }
  sum<--sum
  
  return(sum)}
#### end function


ProbSeq <-function(theta, delta,K){
  
  ### computes probabilities of SEQ model 
  ###  k+1 categories 0,1,  K
  ### delta: k-vector of thresholds
  
  
  probF<-rep(0,K)
  
  ### seq
  for (r in 1:K)probF[r]<-respfctc(theta-delta[r])
  
  ##comp prob
  prob0[1]<-1-probF[1]
  
  Kn1<-K-1
  for (r in 1:Kn1){dum<-1-probF[r+1]  
  for (j in 1:r) dum<-dum*probF[j]
  prob0[r+1]<-dum} #with resp o at first place
  
  dum<-1
  for (j in 1:K) dum<-dum*probF[j]
  prob0[K+1]<-dum
  #sum(prob0)
  
  
  
  return(prob0)}
#### end function



