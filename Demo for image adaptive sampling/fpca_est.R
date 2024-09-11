#Input X is of m*n1*n2*p matrix (number of samples, n1*n2 response function  f(s1,s2), number of sensors)
##d0 is the number of decomposition dimension. in E.q 5
fpca_est<-function(X,sensor_id){
  X1<-X[,,,sensor_id]
  m <- dim(X1)[1]
  n1<- dim(X1)[2]
  n2<- dim(X1)[3]
  #estimate mean-function 
  mua <- apply(X1,2:3,mean)
  #Covariance function in the paper 
  cova = array(0,dim = c(n1,n1,n2,n2))
  for(u in 1:m){
    fd<-array(0,dim = c(n1,n1,n2,n2))
    for(i3 in 1:n1){
      for(i4 in 1:n2){
        fd[,i3,,i4] <- (X1[u,,]-mua)*(X1[u,i3,i4]-mua[i3,i4]) #question
      }
    }
    cova = cova + fd
  }
  cova = cova/m
  cova2 <- array(0, dim = c(n1^2,n2^2))
  for (i in 1:n2){ 
    for (j in 1:n2){
      rowstart<- 1+n2*(i-1)
      rowend <- n2*i
      colstart <- 1+n2*(j-1)
      colend <-n2*j
      cova2[rowstart:rowend,colstart:colend]<-cova[,,i,j]*(1/n1)*(1/n2)
    }
  }
  
  #Eigen Decomposition 
  eio<-eigs_sym(cova2, k = 99, which = "LM")
  eval <- eio$values
  efns <- eio$vectors
  npc <- 0 
  d0 <- 0
  while(npc<0.95){ #specify percentage of covariance 
    d0 <- d0+1
    npc <- sum(eval[1:d0])/sum(eval)
  }
  eval0<-eval[1:d0] #choice of dimension d0 eigenvalue
  efns0<-efns[,1:d0] #choice of dimension d0 eigenvalue
  
  efns_t<-array(0,dim = c(n1,n2,d0))
  for(k in 1:d0){
    efns_tk<- efns0[,k]
    for(j in 1:n2){
      col_s <- 1+n1*(j-1);
      col_d <- n1*j
      efns_t[,j,k] <-  efns_tk[col_s:col_d]
    }
    #normalize eigenfunction to be with norm of 1
    efns_t[,,k]<-efns_t[,,k]*sqrt(sum(efns_t[,,k]*efns_t[,,k]*(1/n1)*(1/n2)))^(-1)
  }
  scores <-array(0,dim = c(m,d0))
  for(i in 1:m){
    scores[i,]<-fpca_score(X1[i,,],efns_t = efns_t,mua = mua)
  }
  covb<-array(0,dim = d0)
  for (k in 1:d0){
    for (i in 1:m){
        demean<- X1[i,,] - mua #100 by 100
        #under selected sensor to fill array
        a<- sum(demean*efns_t[,,k]*(1/n1)*(1/n2))
        b<- sum(demean*efns_t[,,k]*(1/n1)*(1/n2))
        covb[k] <- covb[k] + a*b
    }
    covb[k] <- covb[k]/m
  }
  #output list includes: 
  #1.dimension, 2.template profile, 3.eigenvalue, 4.eigenfunction, 
  #5. mfpca scores #correlation matrix between p sensors. 
  list(d=d0,mua=mua,evalue = eval0,eigenspace = efns_t, scores = scores,covb = covb)
}

#X1 will be a single profile sample readings, so dim(X1) = 1*n1*n2
fpca_score<-function(X,efns_t,mua){
  demeaned <- X-mua #n1*n2
  n1 <- dim(efns_t)[1];n2<-dim(efns_t)[2]
  d <- dim(efns_t)[3]
  scores <-rep(0,length.out = d)
  for(k in 1:d){
      scores[k] = sum(demeaned*efns_t[,,k]*(1/n1)*(1/n2))
  }
  scores #d*1 array
}

