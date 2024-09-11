ARL_initial<-function(testset,obnum,delta,mf_est,q1){ 
  p<-dim(testset)[4]
  #set up initial table and setting up initial compensation
  cusum_score<-rep(0,length.out = p) 
  table <- data.frame(sensor_id = seq(1,p,length.out = p),cusum_score)
  #each score are sum(ksi1+ksi2+...ksi4)
  data_samp <- testset[1,,,]
  ob_init <- m_score(data_samp,mf_est$eigenspace,mf_est$mua) #d*p matrix score table
  for(i in 1:obnum){
    scale <- E1[i,i] #in E.q 10
    table[i,"cusum_score"] <- sum(ob_init[,i])*scale
    table[i,"W_positive"]<-max(umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
   if(obnum <= (p-1)){
   for(i in (obnum+1):p){
     table[i,"cusum_score"] <-NA
     table[i,"W_positive"]<- delta
     table[i,"W_negative"]<- delta
     table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }}
  #order local statistics in Eq. 12
  table["rank"] = rank(-table[,"W_i"],ties.method = "random")
  current_ob <- seq(1,obnum,length.out = obnum)
  next_ob <- subset(table,rank<= obnum)$sensor_id
  #filter the top-q index 
  top_q<- subset(table,rank<= q1)$sensor_id
  list(table = table,current_ob = current_ob,next_ob = next_ob,top_q = top_q)
}

Get_ARL <- function(testset,Corr,L,obnum,delta,mf_est,q1){
  procedure <-vector(mode ="list",length = 500)
  ob_hist<-array(NA,dim = c(500,obnum))
  L_hist <- rep(NA,length.out = dim(testset)[1]) 
  Init <- ARL_initial(testset,obnum,delta,mf_est,q1)
  ob_hist[1,] <- Init$current_ob
  ob_hist[2,] <-Init$next_ob
  procedure[[1]]<-Init$table
  s <- as.matrix(Init$table[,"W_i"])[Init$top_q]
  L_hist[1] <- sqrt(t(s)%*%solve(Corr[Init$top_q,Init$top_q])%*%s)
  if(L_hist[1]>=L){
    return (list(RL = 1,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
  }
  current <- Init$table
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_stage(testset,current,sample,obnum,delta,mf_est,q1)
    next_table<- next_step$table
    s <- as.matrix(next_table[,"W_i"])[next_step$top_q]
    #construct of global statistics in Eq. 13
    L_hist[sample] <- sqrt(t(s)%*%solve(Corr[Init$top_q,Init$top_q])%*%s)
    procedure[[sample]]<-next_table
    if (sample <= (dim(testset)[1]-1)){
      next_ob<-next_step$observed
      ob_hist[sample+1,] <-next_ob
    }
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
    }
    current <- next_table
  }
  list(RL = 500,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist)
}


cusum_stage <- function(Entire,input,example,obnum,delta,mf_est,q1){
  table <- input
  p<-dim(Entire)[4]
  ob_comp <- m_score(Entire[example,,,],mf_est$eigenspace,mf_est$mua) #return d*p
  #random compensation 
  for (i in 1:p){
    if(table[i,]$rank<=obnum){
      #observed sensor and aggregate all the d = 1,2,...d0 for the mfpc scores
      scale <- E1[i,i] #scaling in E.q 10
      table[i,"cusum_score"] = sum(ob_comp[,i])*scale #marginal N(0,1) normal variable
      table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
      table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
      table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
    }
    else{
      table[i,"cusum_score"] <-NA
      table[i,"W_positive"]<-table[i,"W_positive"] + delta
      table[i,"W_negative"]<-table[i,"W_negative"] + delta
      table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
    }
  }
  #order local statistics in Eq. 12
  table$rank <- rank(-1*table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<=obnum)$sensor_id
  #filter the top-q index 
  top_q<- subset(table,rank<=q1)$sensor_id
  list(table = table,observed = ob_sensor,top_q = top_q)
}