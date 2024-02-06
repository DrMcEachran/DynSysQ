#' Rainy River at Manitou Falls: USGS 05133500
#' Daily mean flow from USGS website
#' RAIM from the NWS system, model calibrated based on NWS standard calibration routines: 
#' Anderson, E.A., 2002. Calibration of conceptual hydrologic models for use in river forecasting. Office of Hydrologic Development, US National Weather Service, Silver Spring, MD, 372.
#' 
#' This script will run a long-term simulation (ie non-assimilating of new streamflow data)
#' Training 1950-1999
#' Train linear bias correction based on 2000-2009
#' Test linear bias correction based on 2010-2019

  # read libraries: 
    library(nonlinearTseries)
    library(lubridate)
    library(dplyr)
    library(R2jags)
  
  # knn function: 
    source("Scripts/weightedKNN_longtermForecast.R")
  
  # read in testing and training data: 
    load("Data/manitouTest.RData")
    load("Data/manitouTrain.RData")
  
  # create the Takens vectors: 
    m<-3  # 3 days of streamflow embedding dim
    td<-1 # 1 day time lag
    days.precip<-10  # 7 days of past precip. Will also add one day of perfectly known QPF. 
    days.temp<-3  # 3 days of past observed daily mean temperature.
    
    takSamp<-buildTakens(train$Qcfs, embedding.dim = m, time.lag = td)  # the training attractor for streamflow
    takPred<-buildTakens(test$Qcfs, embedding.dim = m, time.lag = td)  # the testing attractor for streamflow
    
    precipTakTrain<-buildTakens(train$raim,embedding.dim=days.precip,time.lag=1)  # precip training vectors
    precipTakPred<-buildTakens(test$raim, embedding.dim = days.precip, time.lag = 1)  # precip testing vectors
  
    # 1 day of QPF: 
    precipTakTrain<-cbind(precipTakTrain,c(precipTakTrain[2:length(precipTakTrain[,1]),days.precip],10))
    precipTakPred<-cbind(precipTakPred,c(precipTakPred[2:length(precipTakPred[,1]),days.precip],10))
  
    # now the temperature: 
    tempTrain<-buildTakens(train$temp, embedding.dim = days.temp, time.lag = td)
    tempTest<-buildTakens(test$temp, embedding.dim = days.temp, time.lag = td)
    
    # align the vectors temporally: When days.precip is longer than the others. 
    tempTrain<-tempTrain[(days.precip-(days.temp-1)):(length(tempTrain[,1])),] 
    tempTest<-tempTest[(days.precip-(days.temp-1)):(length(tempTest[,1])),] 
    takSamp<-takSamp[(days.precip-(m-1)):(length(takSamp[,1])),] 
    takPred<-takPred[(days.precip-(m-1)):(length(takPred[,1])),] 
  
  # Implement the algorithm: 
    longterm.f<-longterm.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain)
  
  # Fit statistics:
    nse.fcn(test$Qcfs[11:(length(longterm.f)+7)],longterm.f[4:length(longterm.f)]) 
    r2.fcn(test$Qcfs[11:(length(longterm.f)+7)],longterm.f[4:length(longterm.f)]) 
    
    
    
#' Now do a long-term bias correction: 

    # Bias Train here:
    train.bc<-subset(test,year(date)<2014)
    test.bc<-subset(test,year(date)>=2014)  # note: 2014 is the flood of record.
    
    # for the training bc sim: 
      # create the Takens vectors: 

      takSamp<-buildTakens(train$Qcfs, embedding.dim = m, time.lag = td)  # the training attractor for streamflow
      takPred<-buildTakens(train.bc$Qcfs, embedding.dim = m, time.lag = td)  # the testing attractor for streamflow
      
      precipTakTrain<-buildTakens(train$raim,embedding.dim=days.precip,time.lag=1)  # precip training vectors
      precipTakPred<-buildTakens(train.bc$raim, embedding.dim = days.precip, time.lag = 1)  # precip testing vectors
      
      # 1 day of QPF: 
      precipTakTrain<-cbind(precipTakTrain,c(precipTakTrain[2:length(precipTakTrain[,1]),days.precip],10))
      precipTakPred<-cbind(precipTakPred,c(precipTakPred[2:length(precipTakPred[,1]),days.precip],10))
      
      # now the temperature: 
      tempTrain<-buildTakens(train$temp, embedding.dim = days.temp, time.lag = td)
      tempTest<-buildTakens(train.bc$temp, embedding.dim = days.temp, time.lag = td)
      
      # align the vectors temporally: When days.precip is longer than the others. 
      tempTrain<-tempTrain[(days.precip-(days.temp-1)):(length(tempTrain[,1])),] 
      tempTest<-tempTest[(days.precip-(days.temp-1)):(length(tempTest[,1])),] 
      takSamp<-takSamp[(days.precip-(m-1)):(length(takSamp[,1])),] 
      takPred<-takPred[(days.precip-(m-1)):(length(takPred[,1])),] 
      
      train.knnSim<-longterm.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain)
      train.knnSim<-train.knnSim[(m+1):length(train.knnSim)] # forecast starts on day m+1
      train.knnSim<-data.frame("date"=train.bc$date[(days.precip+1):(length(train.knnSim)+days.precip)],"knn"=train.knnSim)
      
    # build the regression: 
      obsvssim<-inner_join(test,train.knnSim)
      load("Data/littleFork.RData")  # little fork river: date, lf=little fork OBSERVED mean daily flows.
      # sim = SIMULATED little fork mean daily flows.
      obsvssim<-inner_join(obsvssim,lf)
      
      
      # Linear Bias Correct: 
      mod.lbc<-function(){
        for(i in 1:n){
          mu[i]<-b0+b1*Qknn[i]+b2*lf[i] 
          Qobs[i]~dnorm(mu[i],prec)
        }
        
        b0~dnorm(0,0.01)
        b1~dnorm(0,0.01)
        b2~dnorm(0,0.01)
        prec<-pow(sd.fdays,-2)
        sd.fdays~dunif(1,20000)
        
        for(j in 1:n){  # simulate values
          Qsims[j]~dnorm(mu[j],prec)
        }
      }
      
      dat.lbc<-list(
        "n"=length(obsvssim$knn),
        "Qobs"=obsvssim$Qcfs,
        "Qknn"=obsvssim$knn,
        "lf"=obsvssim$sim
        
      )
      
      params.lbc<-c("mu","Qsims","b0","b1","b2","sd.fdays")
      
      lbc<-jags.parallel(model.file=mod.lbc, data=dat.lbc, parameters.to.save=params.lbc,
                         n.chains=3, n.iter=4000, n.burnin=1000
      )
      
    # in the true testing time period: 
      takSamp<-buildTakens(train$Qcfs, embedding.dim = m, time.lag = td)  # the training attractor for streamflow
      takPred<-buildTakens(test.bc$Qcfs, embedding.dim = m, time.lag = td)  # the testing attractor for streamflow
      
      precipTakTrain<-buildTakens(train$raim,embedding.dim=days.precip,time.lag=1)  # precip training vectors
      precipTakPred<-buildTakens(test.bc$raim, embedding.dim = days.precip, time.lag = 1)  # precip testing vectors
      
      # 1 day of QPF: 
      precipTakTrain<-cbind(precipTakTrain,c(precipTakTrain[2:length(precipTakTrain[,1]),days.precip],10))
      precipTakPred<-cbind(precipTakPred,c(precipTakPred[2:length(precipTakPred[,1]),days.precip],10))
      
      # now the temperature: 
      tempTrain<-buildTakens(train$temp, embedding.dim = days.temp, time.lag = td)
      tempTest<-buildTakens(test.bc$temp, embedding.dim = days.temp, time.lag = td)
      
      # align the vectors temporally: When days.precip is longer than the others. 
      tempTrain<-tempTrain[(days.precip-(days.temp-1)):(length(tempTrain[,1])),] 
      tempTest<-tempTest[(days.precip-(days.temp-1)):(length(tempTest[,1])),] 
      takSamp<-takSamp[(days.precip-(m-1)):(length(takSamp[,1])),] 
      takPred<-takPred[(days.precip-(m-1)):(length(takPred[,1])),] 
      
      test.knnSim<-longterm.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain)
      test.knnSim<-test.knnSim[(m+1):length(test.knnSim)] # forecast starts on day m+1
      test.knnSim<-data.frame("date"=test.bc$date[(days.precip+1):(length(test.knnSim)+days.precip)],"knn"=test.knnSim)
      load("Data/littleFork.RData")  # little fork river: date, lf=little fork OBSERVED mean daily flows.
      test.knnSim<-inner_join(test.knnSim,lf)
      test.knnSim<-inner_join(test.knnSim,test)
      
    # Now, based on the regression, resample predictions:
      
      sms<-matrix(nrow=5000,ncol=length(test.knnSim$knn))  # sms: for "sims", every row is a timeseries of the bias corrected simulated flows, each row being a sample of the posterior of the parameters
      for(j in 1:5000){
        rand<-ceiling(runif(1,0,3000))
        print(j)
        for(i in 1:length(test.knnSim$knn)){
          mu<-lbc$BUGSoutput$sims.list$b0[rand]+lbc$BUGSoutput$sims.list$b1[rand]*test.knnSim$knn[i]+lbc$BUGSoutput$sims.list$b2[rand]*test.knnSim$sim[i]
          sms[j,i]<-rnorm(1,mu,lbc$BUGSoutput$sims.list$sd.fdays[rand])
        }
      }
      
      # colMeans of sms are the mean estimate
      
      nse.fcn(test.knnSim$Qcfs,colMeans(sms))
      r2.fcn(test.knnSim$Qcfs,colMeans(sms))
      pbias.fcn(test.knnSim$Qcfs,colMeans(sms))
      
      # compared to raw knn:
      nse.fcn(test.knnSim$Qcfs,test.knnSim$knn)
      r2.fcn(test.knnSim$Qcfs,test.knnSim$knn)
      pbias.fcn(test.knnSim$Qcfs,test.knnSim$knn)
      
      bc<-colMeans(sms)
      
# PLOT: FIGURE 8
# ----------
      figure.data<-data.frame("date"=test.knnSim$date,"obs"=test.knnSim$Qcfs,"bc"=bc,"knn"=test.knnSim$knn)
      long.plotThing<-pivot_longer(figure.data,cols=c("obs","bc","knn"),names_to="Timeseries",values_to="Streamflow")
      long.plotThing$Timeseries<-ifelse(long.plotThing$Timeseries=="bc","Bias-Corrected KNN",
                                        ifelse(long.plotThing$Timeseries=="knn","Raw KNN","Observed")
                                        )
      
      rainyBC<-ggplot(data=long.plotThing,aes(date,Streamflow,color=Timeseries))+geom_line(lwd=0.75)+
        theme_bw(base_size=(size=14))+theme(legend.position="right")+
        xlab("Date")+
        ylab("Streamflow [cfs]")+
        scale_color_manual(values=c("red","black","darkgoldenrod4"))
      
      png(filename="RainyRiverStreamflow.png",width=10,height=6,res=300,units="in")
      rainyBC
      dev.off()
      
# ----------
      
      # Run the ENTIRE timeseries in the testing time period 
      # so as to capture long-term trends in the errors using 
      # frequency analysis.

      takSamp<-buildTakens(train$Qcfs, embedding.dim = m, time.lag = td)  # the training attractor for streamflow
      takPred<-buildTakens(test$Qcfs, embedding.dim = m, time.lag = td)  # the testing attractor for streamflow
      precipTakTrain<-buildTakens(train$raim,embedding.dim=days.precip,time.lag=1)  # precip training vectors
      precipTakPred<-buildTakens(test$raim, embedding.dim = days.precip, time.lag = 1)  # precip testing vectors
      precipTakTrain<-cbind(precipTakTrain,c(precipTakTrain[2:length(precipTakTrain[,1]),days.precip],10))
      precipTakPred<-cbind(precipTakPred,c(precipTakPred[2:length(precipTakPred[,1]),days.precip],10))
      tempTrain<-buildTakens(train$temp, embedding.dim = days.temp, time.lag = td)
      tempTest<-buildTakens(test$temp, embedding.dim = days.temp, time.lag = td)
      tempTrain<-tempTrain[(days.precip-(days.temp-1)):(length(tempTrain[,1])),] 
      tempTest<-tempTest[(days.precip-(days.temp-1)):(length(tempTest[,1])),] 
      takSamp<-takSamp[(days.precip-(m-1)):(length(takSamp[,1])),] 
      takPred<-takPred[(days.precip-(m-1)):(length(takPred[,1])),] 
      
      test.longterm.knn<-longterm.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain)
      
      # line up with dates: 
      test.longterm.knn<-test.longterm.knn[(m+1):length(test.longterm.knn)] # forecast starts on day m+1
      test.longterm.knn<-data.frame("date"=test$date[(days.precip+1):(length(test.longterm.knn)+days.precip)],"knn"=test.longterm.knn)
      test.longterm.knn<-inner_join(test.longterm.knn,test)
      
      error<-test.longterm.knn$knn-test.longterm.knn$Qcfs
      perc.error<-(error/test.longterm.knn$Qcfs)*100
      plot(perc.error,type="l")
      abline(h=0,col="red")
      ssp <- spectrum(perc.error,log="no")  
      per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
      plot((1/ssp$freq),ssp$spec,type="l")
      points((1/ssp$freq),ssp$spec)
      abline(v=365/2)
      abline(v=365*4)
      # dominant spikes around 6mos, 4yrs
      
      
      times<-seq(1,length(error))
      reslm <- lm(error ~ sin(2*pi/per*times)+cos(2*pi/per*times))
      summary(reslm)
