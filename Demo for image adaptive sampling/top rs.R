RS_init<-function(testset,obnum,fpca_list,q1){
  p<-dim(testset)[4]
  cusum_score<-rep(0,length=p)
  table <- data.frame(sensor_id = seq(1,p,length.out = p),cusum_score)
  for (i in 1:obnum){
    ob_init<-fpca_score(testset[1,,,i],fpca_list[[i]]$eigenspace,fpca_list[[i]]$mua) 
    table[i,"cusum_score"] <- sum(ob_init)*E[i]
    table[i,"W_positive"]<-max(umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  if(obnum <= (p-1)){
  for(i in (obnum+1):p){
    table[i,"cusum_score"] <-NA
    table[i,"W_positive"]<-0
    table[i,"W_negative"]<-0
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }}
  table["rank"] = rank(-table[,"W_i"],ties.method = "random")
  current_ob <- seq(1,obnum,length.out = obnum)
  next_ob <- sample(1:p,obnum)
  top_q <- subset(table,rank<= q1)$sensor_id
  list(table = table ,current_ob = current_ob,next_ob = next_ob,top_q = top_q)
}


RS_ARL<-function(testset,L,obnum,fpca_list,q1){
  procedure <-vector(mode ="list",length = 500)
  ob_hist<-array(0,dim = c(500,obnum))
  L_hist <- rep(0,dim(testset)[1]) 
  init_list <- RS_init(testset,obnum,fpca_list,q1)
  ob_hist[1,] <- init_list$current_ob
  ob_hist[2,] <- init_list$next_ob
  procedure[[1]]<-init_list$table
  s <- as.matrix(init_list$table[,"W_i"])[init_list$top_q]
  L_hist[1] <- sum(s) 
  if(L_hist[1]>=L){
    return (list(RL = 1,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
  }
  current <- init_list
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_RS(testset,current,sample,obnum,fpca_list,q1) 
    next_table<- next_step$table  
    s <- as.matrix(next_table[,"W_i"])[next_step$top_q]
    L_hist[sample] <- sum(s)
    procedure[[sample]]<-next_table
    if(sample<=(dim(testset)[1]-1)){
      ob_hist[sample+1,] <- next_step$next_ob
    }
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
    }  
    current <- next_step
  }
  list(RL = 500,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist)
}


cusum_RS <- function(Entire,listin,example,obnum,fpca_list,q1){
  p<-dim(Entire)[4]
  table<-listin$table
  for (i in 1:p){
    if(is.element(i,listin$next_ob)){
      ob_comp <- fpca_score(Entire[example,,,i],fpca_list[[i]]$eigenspace,fpca_list[[i]]$mua) #d*1
      table[i,"cusum_score"] <- sum(ob_comp)*E[i] #N(0,1) variable
      table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
      table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
      table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
    }
    else{
      table[i,"cusum_score"] <-NA
      table[i,"W_positive"]<-table[i,"W_positive"]
      table[i,"W_negative"]<-table[i,"W_negative"]
      table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
    }
  }
  table$rank <- rank(-1*table$W_i,ties.method = "random")
  next_ob<- sample(1:p,obnum)
  top_q<-subset(table,rank<= q1)$sensor_id
  list(table = table,next_ob = next_ob,top_q = top_q)
}


