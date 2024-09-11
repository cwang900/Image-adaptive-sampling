#Input X is of m*n1*n2*p matrix (number of samples, n1*n2 response function  f(s1,s2), number of sensors)
##d0 is the number of decomposition dimension. in Eq. 2
mfpca_est<-function(X){
  m <- dim(X)[1]
  n1<- dim(X)[2]
  n2<- dim(X)[3]
  p <- dim(X)[4]
  #estimate mean-function 
  mua <- apply(X,2:4,mean)
  
  #Scalar covariance function in the paper in Eq.A12
  cova = array(0,dim = c(n1,n1,n2,n2))
  for(u in 1:m){
    sum_p = array(0,dim = c(n1,n1,n2,n2))
    for (i in 1:p){
      fd<-array(0,dim = c(n1,n1,n2,n2))
      for(i3 in 1:n1){
        for(i4 in 1:n2){
          fd[,i3,,i4] <- (X[u,,,i]-mua[,,i])*(X[u,i3,i4,i]-mua[i3,i4,i])
        }
      }
      sum_p = sum_p + fd
    }
    cova = cova + sum_p
  }
  cova = cova/m
 
  #arry 100^4 to 10000^2
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
  
  #Eigen Decomepostion 
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
  efns0<-efns[,1:d0] 
  
  #Eigenfunction estimation and reshape in Eqs.A11 and A13
  efns_t<-array(0,dim = c(n1,n2,d0))
  for(k in 1:d0){
    efns_tk<- efns0[,k]
    for(j in 1:n2){
      col_s <- 1+n1*(j-1);
      col_d <- n1*j
      efns_t[,j,k] <-  efns_tk[col_s:col_d]
    }
    #normalize eigenfunction to be with norm of 1 in Eq.A14 of section B of appendix
    efns_t[,,k]<-efns_t[,,k]*sqrt(sum(efns_t[,,k]*efns_t[,,k]*(1/n1)*(1/n2)))^(-1)
  }
  scores <-array(0,dim = c(m,d0,p))
  for(i in 1:m){
    scores[i,,]<-m_score(X[i,,,],efns_t = efns_t,mua = mua)
  }
  #covb in Eq.A16 in appendix B
  covb<-array(0,dim = c(d0,p,p))
  for (k in 1:d0){
    for (i in 1:m){
      a<-rep(NA, length = p)
      b<- rep(NA,length = p)
      for(j in 1:p){
        demean<- X[i,,,] - mua
        #under selected sensor to fill array
        a[j]<- sum(demean[,,j]*efns_t[,,k]*(1/n1)*(1/n2))
        b[j]<- sum(demean[,,j]*efns_t[,,k]*(1/n1)*(1/n2))
      }
      covb[k,,] <- covb[k,,] + a%*%t(b)
    }
  }
  covb <- covb/m
  #output list includes: 
  #1.dimension, 2.template profile, 3.eigenvalue, 4.eigenfunction, 
  #5. mfpca scores 6.correlation matrix between p sensors. 
  list(d=d0,mua=mua,evalue = eval0,eigenspace = efns_t, scores = scores, covb = covb)
}

#X1 will be a single profile sample readings, so dim(X1) = 1*n1*n2*p 
m_score<-function(X1,efns_t,mua){
  demeaned <- X1-mua #n1*n2*p
  n1 <- dim(efns_t)[1];n2<-dim(efns_t)[2]
  p <- dim(mua)[3]
  d <- dim(efns_t)[3]
  #scores
  scores <-array(0,dim = c(d,p))
  #see E.q A15 in appendix B
  for(k in 1:d){
    for(l in 1:p){
    scores[k,l] = sum(demeaned[,,l]*efns_t[,,k]*(1/n1)*(1/n2))
    }
  }
  scores #d*p matrix
}



