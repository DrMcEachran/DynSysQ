#' Outline for assimilating rainy river: 
#' Run the forecasting matrices for the test period
#' The test period will then be split into two periods: one to train the linear bias correction by leading day (2000-2013) and another to test the final results (2014-2019)
#' The bias correction will be according to the knn alone (not Little Fork)
#' 
#' 
#' Note that this script takes a long time to run because at each timestep, we forecast seven days, then the next day we reinitialize and reforecast. 

  # read libraries: 
    library(nonlinearTseries)
    library(lubridate)
    library(dplyr)
    library(R2jags)
  
  # knn function: 
    source("Scripts/weightedKNN_assimilating.R")
  
  # read in testing and training data: 
    load("Data/manitouTest.RData")
    load("Data/manitouTrain.RData")

  # train is the training set, train.bc is the testing period within which we will train the bias correction, and test.bc is where we will test the bias correction.
    
    train.bc<-subset(test,year(date)<2014)
    test.bc<-subset(test,year(date)>=2014)

  # get the prediction matrix. Here, each row is a forecast, with the first m columns being observed data, and the m+1 column being the first forecast ordinate, out to a forecast length of 7 days.   
    
    # create the Takens vectors: 
      m<-3  # 3 days of streamflow embedding dim
      td<-1 # 1 day time lag
      days.precip<-10  # 7 days of past precip. Will also add one day of perfectly known QPF. 
      days.temp<-3  # 3 days of past observed daily mean temperature.
      
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
      
    # Do the train.bc time period to get the bias correction training: 
      pred.matrix.bc<-assimilating.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain,length.forecast=7)
      
    # for the "test.bc" time period:
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
      
      # Do the train.bc time period to get the bias correction training: 
      pred.matrix.test<-assimilating.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain,length.forecast=7)
      
    #' DO THE BIAS CORRECTION USING JAGS: 
    
      # Now that we have the knn matrices, insert the date of the first ordinate. 
      # Remember that the first three columns are observed data, then forecast day one
      # is in column 4, etc. 
    
      # here, the date column represents the first day in the row (observed)
      # the first forecast day is date+days(3)
      preds.bc<-data.frame("date"=train.bc$date[8:(length(pred.matrix.bc[,1])+7)],pred.matrix.bc)
    
      # Now, make this long format for doing a bias correction regression: 
      datfr<-c()   # a LONG format dataframe of the bias-correction time period predictions
      for(i in 1:7){  # one loop thru for each day. 
        current<-data.frame("date"=preds.bc$date+days(i+2),"pred.flow"=preds.bc[,i+4])
        #current<-inner_join(current,lf.sim,by="date")  # was doing littleFork
        #current<-inner_join(current,bigFork,by="date")
        current<-inner_join(current,test,by="date")
        current$fdays<-rep(i,length(current$date))  # fdays is the forecast day (ie days from t0)
        
        datfr<-rbind(datfr,current)  # make the dataframe a LONG format. 
        
      }
      datfr<-data.frame("date"=datfr$date,"fdays"=datfr$fdays,"Qobs"=datfr$Qcfs,"Qknn"=datfr$pred.flow)
      datfr<-na.omit(datfr)
      
      # Create the JAGS regression:
      linreg.biasCorrect<-function(){  
        
        for(j in 1:fdays){
          for(i in (((j-1)*len)+1):(j*len)){
            mu[i]<-b0[j]+b1[j]*Qknn[i]
            Qobs[i]~dnorm(mu[i],prec[j])
          }
        }
        
        for(k in 1:fdays){
          b0[k]~dnorm(0,0.01)  # b0 is the intercept
          b1[k]~dnorm(0,0.01)  # b1 is the slope
          prec[k]<-pow(sd.fdays[k],-2)  # sd.fdays is a standard error for each of the forecast days. Will get bigger on each subsequent day. 
          sd.fdays[k]~dunif(1,50000)
        }
        
        for(j in 1:fdays){  # simulate values according to the posterior. 
          for(i in (((j-1)*len)+1):(j*len)){
            Qsims[i]~dnorm(mu[i],prec[j])
          }
        }
        
      }
      
      dat.linreg<-list(
        "fdays"=7,
        "Qobs"=datfr$Qobs,
        "Qknn"=datfr$Qknn,
        "len"=length(subset(datfr,fdays==1)$date)  # the length of the testing timeseries for each day. Each one is equal.
      )
      
      params<-c("mu","Qsims","b0","b1","sd.fdays")
      
      # run the jags model: 
      linreg.bc.model<-jags.parallel(model.file=linreg.biasCorrect, data=dat.linreg, parameters.to.save=params,
                              n.chains=3, n.iter=4000, n.burnin=1000
      )
      
    #' NOW THAT THE BC REGRESSION IS BUILT, ADJUST THE TEST PERIOD:
      preds.test<-data.frame("date"=test.bc$date[8:(length(pred.matrix.test[,1])+7)],pred.matrix.test)
      preds.test<-na.omit(preds.test)
      
      # Now, make this long format for doing a bias correction regression: 
      datfr.test<-c()   # a LONG format dataframe of the bias-correction time period predictions
      for(i in 1:7){  # one loop thru for each day. 
        current<-data.frame("date"=preds.test$date+days(i+2),"pred.flow"=preds.test[,i+4])
        current<-inner_join(current,test,by="date")
        current$fdays<-rep(i,length(current$date))  # fdays is the forecast day (ie days from t0)
        datfr.test<-rbind(datfr.test,current)  # make the dataframe a LONG format. 
      }
      datfr.test<-data.frame("date"=datfr.test$date,"fdays"=datfr.test$fdays,"Qobs"=datfr.test$Qcfs,"Qknn"=datfr.test$pred.flow)
      datfr.test<-na.omit(datfr.test)
      
      # now, for each day for which we want to assess, create a dataframe. 
      day1<-subset(datfr.test,fdays==1)
      day3<-subset(datfr.test,fdays==3)
      day7<-subset(datfr.test,fdays==7)
    
    numSims=10000
    warmState.estimate.day1<-matrix(nrow=length(day1$fdays),ncol=numSims)   # columns are each posterior sample, and rows are each day. 
    warmState.estimate.day3<-matrix(nrow=length(day3$fdays),ncol=numSims)
    warmState.estimate.day7<-matrix(nrow=length(day7$fdays),ncol=numSims)
    
    for(j in 1:numSims){
      print(j)
      rand<-ceiling(runif(1,0,3000))  # 3000 McMC samples are saved off from the posterior, so get a random number between 1 and 3000
      for(i in 1:length(day1$fdays)){  # each forecast day matrix is equal length
        
        mu.day1<-linreg.bc.model$BUGSoutput$sims.list$b0[rand,1]+linreg.bc.model$BUGSoutput$sims.list$b1[rand,1]*day1$Qknn[i]
        warmState.estimate.day1[i,j]<-rnorm(1,mu.day1,linreg.bc.model$BUGSoutput$sims.list$sd.fdays[rand,1])
        
        mu.day3<-linreg.bc.model$BUGSoutput$sims.list$b0[rand,3]+linreg.bc.model$BUGSoutput$sims.list$b1[rand,3]*day3$Qknn[i]
        warmState.estimate.day3[i,j]<-rnorm(1,mu.day3,linreg.bc.model$BUGSoutput$sims.list$sd.fdays[rand,3])
        
        mu.day7<-linreg.bc.model$BUGSoutput$sims.list$b0[rand,7]+linreg.bc.model$BUGSoutput$sims.list$b1[rand,7]*day7$Qknn[i]
        warmState.estimate.day7[i,j]<-rnorm(1,mu.day7,linreg.bc.model$BUGSoutput$sims.list$sd.fdays[rand,7])
        
      }
    }
    
    # get the uncertainty:
    
    for(i in 1:length(warmState.estimate.day7[,1])){
      # note: have to use as.numeric to get rid of column quantile label. 
      day1$up[i]<-quantile(warmState.estimate.day1[i,],0.975) # store the 95% credible interval upper bounds here
      day1$lo[i]<-quantile(warmState.estimate.day1[i,],0.025) # store the lower bounds here
      day1$bc.mean[i]<-mean(warmState.estimate.day1[i,])  # bias corrected mean value is bc.mean
      
      day3$up[i]<-quantile(warmState.estimate.day3[i,],0.975)
      day3$lo[i]<-quantile(warmState.estimate.day3[i,],0.025)
      day3$bc.mean[i]<-mean(warmState.estimate.day3[i,])
      
      day7$up[i]<-quantile(warmState.estimate.day7[i,],0.975)
      day7$lo[i]<-quantile(warmState.estimate.day7[i,],0.025)
      day7$bc.mean[i]<-mean(warmState.estimate.day7[i,])

    }
    
  # summary stats for warm state estimates, no bias correction: 
    nse.fcn(day1$Qobs,day1$Qknn)
    nse.fcn(day3$Qobs,day3$Qknn)
    nse.fcn(day7$Qobs,day7$Qknn)
    
    r2.fcn(day1$Qobs,day1$Qknn)
    r2.fcn(day3$Qobs,day3$Qknn)
    r2.fcn(day7$Qobs,day7$Qknn)
    
    pbias.fcn(day1$Qobs,day1$Qknn)
    pbias.fcn(day3$Qobs,day3$Qknn)
    pbias.fcn(day7$Qobs,day7$Qknn)
    
  # summary stats for warm state estimates, WITH bias correction: 
    nse.fcn(day1$Qobs,day1$bc.mean)
    nse.fcn(day3$Qobs,day3$bc.mean)
    nse.fcn(day7$Qobs,day7$bc.mean)
    
    r2.fcn(day1$Qobs,day1$bc.mean)
    r2.fcn(day3$Qobs,day3$bc.mean)
    r2.fcn(day7$Qobs,day7$bc.mean)
    
    pbias.fcn(day1$Qobs,day1$bc.mean)
    pbias.fcn(day3$Qobs,day3$bc.mean)
    pbias.fcn(day7$Qobs,day7$bc.mean)    
  
  
  #' ASSESSMENT: FLOOD-OF-RECORD: 
  
    # PLOT ONE: PROBABILITY DENSITIES:
    png(filename="RainyProbabilities.png",width=10,height=6,res=300,units="in")
    plot(density(warmState.estimate.day1[which.max(day1$Qobs),]),
         xlab="Streamflow [cfs]",main="",
         ylab="Density",
         xlim=c(40000,82000),lwd=4,cex.axis=2,cex.lab=1.75
    )
    abline(v=day1$Qobs[which.max(day1$Qobs)],lty=2,col="red",lwd=2)
    lines(density(warmState.estimate.day7[which.max(day7$Qobs),]),
          col="darkgoldenrod4",lwd=4
    )
    lines(density(warmState.estimate.day3[which.max(day3$Qobs),]),
          col="blue",lwd=4
    )
    legend("topleft",c("7-day lead time","3-day lead time","1-day lead time","Actual Peak Value"),lty=c(1,1,1,1,2),
           lwd=c(2,2,2),col=c("darkgoldenrod4","blue","black","red"),cex=1.5
    )
    dev.off()
    
    # PLOT TWO: THE "OPERATIONAL" FORECAST PICTURE LEADING UP TO THE FLOOD-OF-RECORD: 
    # Hand-pick out the area around the flood-of-record
    # remember, back to preds.test: 
    # X4 is day1
    # X6 is day3
    # X10 is day7

    # REMEMBER X3=t0
    # but, the date column has the date of t-2 ... So to use this to get the date, need the location of the
    # flood of record in the X1 column
    loc.max<-which.max(preds.test$X1)
    
    d1.pred<-preds.test[preds.test$date==preds.test[loc.max,]$date-days(3),]  # the prediction with one day lead time before the flood of record
    d3.pred<-preds.test[preds.test$date==preds.test[loc.max,]$date-days(5),]  # the prediction with 3 days lead time, etc. 
    d7.pred<-preds.test[preds.test$date==preds.test[loc.max,]$date-days(9),]
    d0.pred<-preds.test[preds.test$date==preds.test[loc.max,]$date-days(2),] # the forecast of the day-of the flood of record
    

    
    # for the plot, I want to show the intialization day that is observed blend up into the kNN forecast. 
    # the dates will range from date_X3 to date_X3+days(7) 
    
    pk2014<-subset(test,year(date)==2014)
    
    png(filename="RainyTimeseries.png",width=10,height=6,res=300,units="in")
    
    plot(pk2014$date[157:178],pk2014$Qcfs[157:178],col="red",xlab="date",ylab="Flow [cfs]",cex.axis=1.75,cex.lab=1.75,lwd=2,ylim=c(40000,75000))
    lines(pk2014$date[125:200],pk2014$Qcfs[125:200],lwd=2,col="red")
    
    points(seq.Date(from=preds.test$date[loc.max-1],to=preds.test$date[loc.max+6],by="day"),d1.pred[,4:11],col="black",pch=16)
    lines(seq.Date(from=preds.test$date[loc.max-1],to=preds.test$date[loc.max+6],by="day"),d1.pred[,4:11],col="black",lwd=4)
    
    points(seq.Date(from=preds.test$date[loc.max-3],to=preds.test$date[loc.max+4],by="day"),d3.pred[,4:11],col="blue",pch=16)
    lines(seq.Date(from=preds.test$date[loc.max-3],to=preds.test$date[loc.max+4],by="day"),d3.pred[,4:11],col="blue",lwd=4)
    
    lines(seq.Date(from=preds.test$date[loc.max-7],to=preds.test$date[loc.max],by="day"),d7.pred[,4:11],col="darkgoldenrod4",lwd=4)
    points(seq.Date(from=preds.test$date[loc.max-7],to=preds.test$date[loc.max],by="day"),d7.pred[,4:11],col="darkgoldenrod4",pch=16)
    
    points(seq.Date(from=preds.test$date[loc.max],to=preds.test$date[loc.max+7],by="day"),d0.pred[,4:11],col="darkgray",pch=16)
    lines(seq.Date(from=preds.test$date[loc.max],to=preds.test$date[loc.max+7],by="day"),d0.pred[,4:11],col="darkgray",lwd=4)
    
    abline(v=as.Date("2014-06-17"),lwd=2,lty=2,col="red")
    
    legend("topleft",c("7-day lead time","3-day lead time","1-day lead time","Actual Peak Value","0-day:Recession"),lty=c(1,1,1,2,1),
           lwd=c(2,2,2,2,2),col=c("darkgoldenrod4","blue","black","red","darkgray"),cex=1.5
    )
    
    dev.off()
    
    
    