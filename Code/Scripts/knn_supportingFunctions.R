
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


pw.knn<-function(k,par,pp,nn.pp){  # Sum of squares to minimize
  # k is the k nearest neighbors (rows = neighbor vectors, cols = m+1 (the m past readings and the m+1 point it is traveling to))
  # par are the parameters for the trajectories
  # pp is the current precip vector
  # pp.nn is the precip vector associated with each nearest neighbor
  alpha<-par[1]
  betas<-par[2:length(par)]
  
  d<-c()
  er<-c()
  w<-c()
  minimizeMe<-c()
  
  # first find the distances: 
  for(i in 1:length(k[,1])){
    d[i]<-euc.dist(pp,nn.pp[i,])
  }
  
  if(sum(d)==0){
    d<-rep(1,length(d))
  }
  
  d=d+1  # add some arbitraty value so that we don't blow up at distance = 0 for the weight
  
  b=mean(d) # some value
  w=exp(-0.5*((d/b)^2))  # Gaussian kernel weighting
  #w=w/sum(w)
  #w = ifelse(d<b,(1-((d/b)^2))^2,0)
  
  
  #d<-(1-(d/max(d)))
  
  for(i in 1:length(k[,1])){
    
    er[i]<-(w[i])*(((alpha+(betas%*%k[i,1:length(k[1,])-1]))-k[i,length(k[1,])])^2.0)
    minimizeMe[i]<-er[i]
  }
  min<-sum(minimizeMe)
  return(min)
}

knn<-function(k,par){  # Sum of squares to minimize
  alpha<-par[1]
  betas<-par[2:length(par)]
  minimizeMe<-c()
  for(i in 1:length(k[,1])){
    minimizeMe[i]<-(((alpha+(betas%*%k[i,1:(length(k[1,])-1)]))-k[i,length(k[1,])])^2.0)
  }
  mum<-sum(minimizeMe)
  return(mum)
}



# GOF functions:
nse.fcn<-function(obs,sim){
  numerator<-sum((obs-sim)^2)
  denominator<-sum((obs-mean(obs))^2)
  calc<-(1-(numerator/denominator))
  return(calc)
}

r2.fcn<-function(obs,sim){
  calc<-((cor(obs,sim))^2)
  return(calc)
}

rsr.fcn<-function(obs,sim){
  calc<-sqrt(sum((obs-sim)^2))/sqrt(sum((obs-mean(obs))^2))
  return(calc)
}

pbias.fcn<-function(obs,sim){
  calc<-sum((sim-obs)*100)/sum(obs)
  return(calc)
}

