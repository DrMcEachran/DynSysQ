assimilating.forecaster<-function(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain,length.forecast=7){
  source("Scripts/knn_supportingFunctions.R")

  # pred.matrix: 
  # each row corresponds to initializing on a Takens vector. Each column corresponds to a time ordinate. 
  # the first m columns are the starting vector (each column one of the dimensions)
  # Starting at m+1, each column is a future time vector. column m is t0. Thru a total length of length.forecast+1
  
  pred.matrix<-matrix(nrow=length(takPred[,1]),ncol=(m+length.forecast))
  
  for(q in 1:(length(takPred[,1])-(length.forecast+1))){
    
    forecast<-c(takPred[q,]) # initialize the forecast with the takens vector
    
    pp<-c(    # day q through q+days+1 (ie one day qpf), then after 1 day qpf go with just snowmelt. 
      precipTakPred[q,],test$snowmelt[(q+days.precip+1):(q+days.precip+length.forecast)]
    )
    pp<-buildTakens(pp,days.precip+1,1) # from the vector: c(observed precip, 1-day qpf, and then snowmelt in the forecast period)
    tt<-tempTest[q:(q+length.forecast),] # for building the forecast temps: perfectly known
    
    print(q)
    
    for(j in m:(m+(length.forecast-1))){
      dist<-c()
      precip.dist<-c()
      temp.dist<-c()
      for(i in 1:length(takSamp[,1])){
        dist[i]<-euc.dist(forecast[(j-(m-1)):j],takSamp[i,])
        precip.dist[i]<-euc.dist(pp[j-m+1,],precipTakTrain[i,])
        temp.dist[i]<-euc.dist(tt[j-m+1,],tempTrain[i,])
      }
      
      dist<-(dist-mean(dist))/sd(dist)
      precip.dist<-(precip.dist-mean(precip.dist))/sd(precip.dist)
      temp.dist<-(temp.dist-mean(temp.dist))/sd(temp.dist)
      dist<-dist+precip.dist+temp.dist
      
      kSamp<-matrix(nrow=2*m,ncol=m+1)
      pSamp<-matrix(nrow=2*m,ncol=days.precip+1)  # one day of forecast precipitation/snowmelt
      
      locs<-c()
      for(i in 1:length(kSamp[,1])){  # fill in the k nearest neighbors
        locs[i]<-which.min(dist[1:(length(dist)-1)])
        kSamp[i,] <- c(takSamp[locs[i],],takSamp[locs[i]+1,m])
        dist[locs[i]]<-max(dist)+1.0
        pSamp[i,] <- precipTakTrain[locs[i],]
      }

      lin.reg<-optim(par=c(1,rep(1,m)),fn=pw.knn,  method="L-BFGS-B",lower=0, hessian=F, k=kSamp,pp=pp[j-m+1],nn.pp=pSamp)
      
      pred.inits<-forecast[(j-(m-1)):j]
      al<-lin.reg$par[1]
      be<-lin.reg$par[2:length(lin.reg$par)]   
      pred<-al+be%*%pred.inits
      forecast[j+1]<-pred
    }
    pred.matrix[q,]<-forecast
  }

  return(pred.matrix)
  
}
