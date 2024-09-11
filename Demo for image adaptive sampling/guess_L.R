library(MASS)

mfpca_L<-function(L,obnum,model,delta,q1){
  oc_run<-array(NA,dim = 500)
  for(i in 1:500){
    x_ic<-Gen_X(500)$data
    oc_run[i]<-Get_ARL(x_ic,Corr,L,obnum,delta = delta,model,q1)$RL
  }
  mean(oc_run)
}

fpca_L<-function(L,obnum,fpca_list,delta,q1){
  oc_run<-array(NA,dim = 500)
  for(i in 1:500){
    x_ic<-Gen_X(500)$data
    oc_run[i]<-Get_ARLi(x_ic,L,obnum,delta,fpca_list,q1)$RL
  }
  mean(oc_run)
}

RS_L<-function(L,obnum,fpca_list,q1){
  oc_run<-array(NA,dim = 500)
  for(i in 1:500){
    x_ic<-Gen_X(500)$data
    oc_run[i]<-RS_ARL(x_ic,L,obnum,fpca_list,q1)$RL
  }
  mean(oc_run)
}

mfpca_train<-function(model,q1){
  newL<-array(NA,dim = c(3,3))
  Rlist<-c(10,7,4)
  delta_list<-c(0.01,0.05,0.1)
  #train for MFPCA
  for(i in seq(Rlist)){
    for(j in seq(delta_list)){
      high <-100; low<-0; curr_L<-50; 
      ARL0<-mfpca_L(curr_L,obnum = Rlist[i],model,delta = delta_list[j],q1)
      while(abs(ARL0-200)>=5){#ensure ARL0 is close to 200 
        if(ARL0>= 205){
          high <- curr_L
          curr_L<-(high+low)/2
          ARL0<-mfpca_L(curr_L,obnum = Rlist[i],model,delta = delta_list[j],q1)
        }
        else if(ARL0<= 195){
          low <- curr_L
          curr_L<-(high+low)/2
          ARL0<-mfpca_L(curr_L,obnum = Rlist[i],model,delta = delta_list[j],q1)
        }
      }
      newL[i,j]<-curr_L
      rm(ARL0,high,low,curr_L)
    }
  }
  return(newL)
}

full_train<-function(model,q1){
  #train for MFPCA full
   high <-100; low<-0; curr_L<-50; 
   ARL0<-mfpca_L(curr_L,obnum = 25,model,delta = NA ,q1)
   while(abs(ARL0-200)>=5){#ensure ARL0 is close to 200 
        if(ARL0>= 205){
          high <- curr_L
          curr_L<-(high+low)/2
          ARL0<-mfpca_L(curr_L,obnum = Rlist[i],model,delta = delta_list[j],q1)
        }
        else if(ARL0<= 195){
          low <- curr_L
          curr_L<-(high+low)/2
          ARL0<-mfpca_L(curr_L,obnum = Rlist[i],model,delta = delta_list[j],q1)
          }}
  newL<-curr_L
  return(newL)
}

fpca_train<-function(fpca_list,q1){
  newL<-array(NA,dim = c(3,3))
  Rlist<-c(10,7,4)
  delta_list<-c(0.01,0.05,0.1)
  #train for FPCA
  for(i in seq(Rlist)){
    for(j in seq(delta_list)){
      high <-100; low<-0; curr_L<-50;
      ARL0<-fpca_L(curr_L,obnum = Rlist[i],fpca_list,delta = delta_list[j],q1)
      while(abs(ARL0-200)>=5){#ensure ARL0 is close to 200 
        if(ARL0>= 205){
          high <- curr_L
          curr_L<-(high+low)/2
          ARL0<-fpca_L(curr_L,obnum = Rlist[i],fpca_list,delta = delta_list[j],q1)
        }
        else if(ARL0<= 195){
          low <- curr_L
          curr_L<-(high+low)/2
          ARL0<-fpca_L(curr_L,obnum = Rlist[i],fpca_list,delta = delta_list[j],q1)
        }
      }
      newL[i,j]<-curr_L
      rm(ARL0,high,low,curr_L)
    }
  }
  return(newL)
}

RS_train<-function(Ls,fpca_list,q1){
  newL<-rep(NA,3)
  Rlist<-c(10,7,4)
  #train for RS
  for(i in seq(Rlist)){
    high <-100; low<-0; curr_L<-50; 
    ARL0<-RS_L(curr_L,obnum = Rlist[i],fpca_list,q1)
    while(abs(ARL0-200)>=5){#ensure ARL0 is close to 200 
      if(ARL0>= 205){
        high <- curr_L
        curr_L<-(high+low)/2
        ARL0<-RS_L(curr_L,obnum = Rlist[i],fpca_list,q1)
      }
      else if(ARL0<= 195){
        low <- curr_L
        curr_L<-(high+low)/2
        ARL0<-RS_L(curr_L,obnum = Rlist[i],fpca_list,q1)
      }
    }
    newL[i]<-curr_L
    rm(ARL0,high,low,curr_L)
  }
  #return newL
  return(newL)
}


#Example on estimating control limits h
Ls_model1<-array(NA,dim = c(3,7))
m_Ls<-mfpca_train(model,q1 = 4); Ls_model1[,1:3]<-m_Ls
Lsfull1<-full_train(model,q1 = 4)
f_Ls<-fpca_train(fpca_list,q1 = 4); Ls_model1[,4:6]<-f_Ls
R_Ls<-RS_train(fpca_list,q1 = 4); Ls_model1[,7]<-R_Ls

