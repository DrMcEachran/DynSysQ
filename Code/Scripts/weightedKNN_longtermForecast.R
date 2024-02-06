# create the longterm forecast function. 

longterm.forecaster<-function(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain){
  source("Scripts/knn_supportingFunctions.R")
  
  m<-length(takSamp[1,])
  days.precip<-length(precipTakPred[1,])
  
  longterm.f<-c(takPred[1,])
  
  for(j in m:(length(takPred[,1]))){    #length(takPred[,1])){
    print(j)
    dist<-c()
    precip.dist<-c()
    temp.dist<-c()
    for(i in 1:length(takSamp[,1])){
      dist[i]<-euc.dist(longterm.f[(j-(m-1)):j],takSamp[i,])
      precip.dist[i]<-euc.dist(precipTakPred[(j-m+1),],precipTakTrain[i,])
      temp.dist[i]<-euc.dist(tempTest[(j-m+1),],tempTrain[i,])
    }
    
    dist<-(dist-mean(dist))/sd(dist)
    precip.dist<-(precip.dist-mean(precip.dist))/sd(precip.dist)
    temp.dist<-(temp.dist-mean(temp.dist))/sd(temp.dist)
    dist<-dist+precip.dist+temp.dist
    
    kSamp<-matrix(nrow=2*m,ncol=m+1)
    pSamp<-matrix(nrow=2*m,ncol=days.precip)  # ** without QPF: matrix(nrow=2*m,ncol=m)
    locs<-c()
    for(i in 1:length(kSamp[,1])){  # fill in the k nearest neighbors
      locs[i]<-which.min(dist[1:(length(dist)-1)])
      kSamp[i,] <- c(takSamp[locs[i],],takSamp[locs[i]+1,m])
      dist[locs[i]]<-max(dist)+1.0
      pSamp[i,] <- precipTakTrain[locs[i],]
    }
    
    lin.reg<-optim(par=c(1,rep(1,m)),fn=pw.knn,  method="L-BFGS-B",lower=0, hessian=F, k=kSamp,pp=precipTakPred[(j-m+1),],nn.pp=pSamp)
    
    # Now that our regression is defined, we want to predict the next value:
    pred.inits<-longterm.f[(j-(m-1)):j]
    al<-lin.reg$par[1]     
    be<-lin.reg$par[2:length(lin.reg$par)]       
    pred<-al+be%*%pred.inits
    longterm.f[j+1]<-pred
  }
  return(longterm.f)
}

