library(MASS)
library(RSpectra)

source('mfpca_estimation.R')
source('fpca_est.R')
source('top arl.R')
source('top fpca.R')
source('top rs.R')

Gen_X<-function(m){
  l=seq(1,25,length.out = 25)
  n=length(l)
  Sigma=matrix(data=NA, nrow=n, ncol=n)
  for(i in 1:n){
    for (j in 1:n){
      Sigma[i,j]= 0.8^abs(i-j) #model II signal
    }
  }
  s1<-seq(0.02,1,length.out = 50)
  s2<-seq(0.02,1,length.out = 50)
  vk <-array(NA,dim = c(50,50,4))
  
  vk[,,1]<-2*cos(2*pi*s1)%*%t(cos(2*pi*s2))
  vk[,,2]<-2*sin(2*pi*s1)%*%t(sin(2*pi*s2))
  vk[,,3]<-2*cos(4*pi*s1)%*%t(cos(4*pi*s2))
  vk[,,4]<-2*sin(4*pi*s1)%*%t(sin(4*pi*s2))
  
  X <-array(NA,dim = c(m,50,50,25))
  for (i in 1:m){
    score_k <- mvrnorm(4,mu = rep(0,25),Sigma = Sigma)
    part<-array(rep(0,length = 50*50*25),dim = c(50,50,25))
    for(d in 1:4){
      part1<-array(rep(0,length = 50*50*25), dim = c(50,50,25))
      for(j in 1:25){
        part1[,,j]<-vk[,,d]*score_k[d,j] #add mean function 
      }
      part = part + part1 
    }
    X[i,,,] = part
  }
  list(data = X,basis = vk) #basis is also attached 
}

Gen_OC<-function(m,mu_degree,oc_channel){
  l=seq(1,25,length.out = 25)
  n=length(l)
  Sigma=matrix(data=NA, nrow=n, ncol=n)
  for(i in 1:n){
    for (j in 1:n){
      Sigma[i,j]= 0.8^abs(i-j) #model II signal
    }
  }
  s1<-seq(0.02,1,length.out = 50)
  s2<-seq(0.02,1,length.out = 50)
  vk <-array(NA,dim = c(50,50,4))
  vk[,,1]<-2*cos(2*pi*s1)%*%t(cos(2*pi*s2))
  vk[,,2]<-2*sin(2*pi*s1)%*%t(sin(2*pi*s2))
  vk[,,3]<-2*cos(4*pi*s1)%*%t(cos(4*pi*s2))
  vk[,,4]<-2*sin(4*pi*s1)%*%t(sin(4*pi*s2))
  
  X <-array(NA,dim = c(m,50,50,25))
  for (i in 1:m){
    score_k <- mvrnorm(4,mu = rep(0,25),Sigma = Sigma)
    part<-array(0,dim = c(50,50,25))
    for(d in 1:4){
      part1<-array(0, dim = c(50,50,25))
      for(j in 1:25){
        part1[,,j]<-vk[,,d]*score_k[d,j]
      }
      part = part + part1 
    }
    #adjust the mean-function
    u<-array(0,dim = c(50,50,25))
    shiftfun<-2*exp(s1+0.5)%*%t(exp(s2+0.5)) #out-of-control setting I
    for (j in oc_channel){
      u[,,j] <-mu_degree*shiftfun
    }
    X[i,,,] = part +  u
  }
  X #basis is also attached 
}

start<-Sys.time()

#training models for 2D MFPCA and FPCA
set.seed(502)
X_param<-Gen_X(200)$data
model<-mfpca_est(X_param) #exp_dim4
fpca_list<-vector(mode = 'list',length = 25)
for(i in 1:25){
  fpca_list[[i]]<-fpca_est(X_param,sensor_id = i)
}

#estimate Corr and Scale the Scores in Eq.10
p<-dim(X_param)[4]
E <-rep(0,length = p)
for(j in 1:p){
  Covb = sum(fpca_list[[j]]$covb)
  E[j] <- 1/sqrt(Covb)
}
rm(Covb)
Covb_mfpca <- array(0,dim = c(p,p)) #aggregate all the d = 1,2,...d0 for the covariance matrix 
for (i in 1:model$d){
  Covb_mfpca <- Covb_mfpca + model$covb[i,,]
}
E1 <- solve(diag(sqrt(diag(Covb_mfpca)))) #Estimated correlation matrix after scaling
Corr <- E1%*%Covb_mfpca%*%E1
rm(p,Covb_mfpca,i)
end<-Sys.time();end-start;

#control limits estimated in "guess_L.R" with corresponded data generated from model II
Ls_model2<-array(NA,dim = c(3,7))
#Column 1-3 MFPCA delta = 0.01,0.05,0.1, Column 4-6 FPCA delta = 0.01,0.05,0.1 and column 7 is FPCA RS
Ls_model2[1,]<-c(8.971524,9.046043,9.271012,15.82838,16.31462,16.73709,17.50471)
Ls_model2[2,]<-c(8.737776,9.011828,9.166346,14.01775,14.80486,16.08396,17.15413)
Ls_model2[3,]<-c(8.578908,8.762241,9.161847,11.01116,12.88190,16.53709,15.99127)
Lsfull2 <-9.34 #control limits for MFPCA fully observed

new_fillans<-function(Ls,obnum,fpca_list,model,q1){
  del <- c(0.01,0.05,0.1)
  oc_run<-array(0,dim = c(6,7,500)) #rep = 500
  oc_sensor<-c(9,16,17,18) 
  #oc_sensor<-c(16,17)
  mu_degree<-c(0.8,1.2,1.6,2.0,3.0,4.0)
  for(i in 1:500){
    for(m in 1:6){
      x_oc<-Gen_OC(500,mu_degree = mu_degree[m],oc_channel = oc_sensor)
      for(j in 1:3){
        oc_run[m,j,i]<-Get_ARL(x_oc,Corr,Ls[j],obnum,delta = del[j],model,q1)$RL
        oc_run[m,j+3,i]<-Get_ARLi(x_oc,Ls[j+3],obnum,delta = del[j],fpca_list,q1)$RL
      }
      oc_run[m,7,i]<-RS_ARL(x_oc,Ls[7],obnum,fpca_list,q1)$RL
    }
  }
  ans<-array(NA,dim = c(12,7))
  for(m in 1:6){
    ans[(2*m-1),]<-apply(oc_run[m,,],1,mean)
    ans[2*m,]<-apply(oc_run[m,,],1,sd)/sqrt(500)
  }
  ans
}

new_full<-function(Ls,obnum,fpca_list,model,q1){
  oc_run<-array(0,dim = c(7,500)) #rep = 500
  oc_sensor<-c(9,16,17,18)
  mu_degree<-c(0,0.8,1.2,1.6,2.0,3.0,4.0) #exp-shift
  for(i in 1:500){ #500
    for(m in 1:7){
      x_oc<-Gen_OC(500,mu_degree = mu_degree[m],oc_channel = oc_sensor)
      #fill oc_run
      oc_run[m,i]<-Get_ARL(x_oc,Corr,Ls[1],obnum,delta = NA,model,q1)$RL
    }
  }
  ans<-array(NA,dim = c(14,1))
  for(m in 1:7){
    ans[(2*m-1),]<-mean(oc_run[m,])
    ans[2*m,]<-sd(oc_run[m,])/sqrt(500) #adjust based on rep = 500
  }
  ans
}

umin<-1.5
output2<-array(NA,dim = c(36,7))
output2[1:12,]<-new_fillans(Ls = Ls_model2[1,],obnum = 10,fpca_list,model, q1 = 4)
output2[13:24,]<-new_fillans(Ls = Ls_model2[2,],obnum = 7,fpca_list,model, q1 = 4)
output2[25:36,]<-new_fillans(Ls = Ls_model2[3,],obnum = 4,fpca_list,model, q1 = 4)
full2<-new_full(Ls = Lsfull2,obnum = 25,fpca_list,model, q1 = 4)

