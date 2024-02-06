#' COWEETA HYDROLOGIC LAB
#' Using the small reference catchment WS18, and long-term rain gage 06 near the 
#' main "visitor center" area. Not on the watershed but nearby.
#' RRG06 (1936-2019)
#' WS18 (1937-2021)
#'     12.5ha
#'     

  # read libraries: 
    library(nonlinearTseries)
    library(lubridate)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  
  # knn function: 
    source("Scripts/weightedKNN_longtermForecast.R")

  # read in testing and training data: 
    load("Data/CoweetaTest.RData")
    load("Data/CoweetaTrain.RData")
  
  # create the Takens vectors: 
    m<-3  # 3 days of streamflow embedding dim
    td<-1 # 1 day time lag
    days.precip<-7  # 7 days of past precip. Will also add one day of perfectly known QPF. 
    days.temp<-3  # 3 days of past observed daily mean temperature.
    
    takSamp<-buildTakens(train$Qmm, embedding.dim = m, time.lag = td)  # the training attractor for streamflow
    takPred<-buildTakens(test$Qmm, embedding.dim = m, time.lag = td)  # the testing attractor for streamflow
    
    precipTakTrain<-buildTakens(train$Pmm,embedding.dim=days.precip,time.lag=1)  # precip training vectors
    precipTakPred<-buildTakens(test$Pmm, embedding.dim = days.precip, time.lag = 1)  # precip testing vectors
    
    # 1 day of QPF: 
    precipTakTrain<-cbind(precipTakTrain,c(precipTakTrain[2:length(precipTakTrain[,1]),days.precip],10))
    precipTakPred<-cbind(precipTakPred,c(precipTakPred[2:length(precipTakPred[,1]),days.precip],10))
    
    # now the temperature: 
    tempTrain<-buildTakens(train$tmean, embedding.dim = days.temp, time.lag = td)
    tempTest<-buildTakens(test$tmean, embedding.dim = days.temp, time.lag = td)
  
    # align the vectors temporally: When days.precip is longer than the others. 
    tempTrain<-tempTrain[(days.precip-(days.temp-1)):(length(tempTrain[,1])),] 
    tempTest<-tempTest[(days.precip-(days.temp-1)):(length(tempTest[,1])),] 
    takSamp<-takSamp[(days.precip-(m-1)):(length(takSamp[,1])),] 
    takPred<-takPred[(days.precip-(m-1)):(length(takPred[,1])),] 
    
  # Implement the algorithm: 
    longterm.f<-longterm.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain)

  # Fit statistics:
    nse.fcn(test$Qmm[8:(length(longterm.f)+4)],longterm.f[4:length(longterm.f)]) 
    r2.fcn(test$Qmm[8:(length(longterm.f)+4)],longterm.f[4:length(longterm.f)])
    
  # Plot streamflow data: 
    coweeta<-data.frame("date"=test$date[8:(length(longterm.f)+4)],"obs"=test$Qmm[8:(length(longterm.f)+4)],"sim"=longterm.f[4:length(longterm.f)])
    
    c.0506<-subset(coweeta,year(date)>=2005 & year(date)<=2006)   # only plot 2005 and 2006 just to get a closer look and better plot. 
    
    long.forPlot<-pivot_longer(c.0506,cols=c("obs","sim"),names_to="Timeseries",values_to="Streamflow")
    
    cow05<-ggplot(data=long.forPlot,aes(date,Streamflow,color=Timeseries))+geom_line(lwd=0.75)+
      theme_bw(base_size=(size=20))+
      #ylim(0,20)+
      xlab("Date")+
      ylab("Streamflow [mm/day]")+
      scale_color_manual(values=c("black","red"))
    cow05
    
    png(filename="coweeta_SF.png",width=10,height=6,res=300,units="in")
    cow05
    dev.off()
    
  # ERROR ANALYSIS: 
    coweeta$er.percent<-((coweeta$sim-coweeta$obs)/coweeta$obs)*100  # add percent error
    er.plot<-ggplot(data=coweeta,aes(date,er.percent))+geom_line(lwd=0.75)+
      theme_bw(base_size=(size=20))+
      #ylim(0,20)+
      xlab("Date")+
      ylab("Model Error [%]")+
      scale_color_manual(values=c("black"))+
      geom_abline(intercept=0,slope=0,color="red",lwd=1.5)
    er.plot
    
    # some percentage errors = Infinity (obs = 0), so replace with NA as these are rare so won't really disrupt spectral analysis
    coweeta$er.percent<-ifelse(coweeta$er.percent==Inf, 0,coweeta$er.percent)
    
    ssp <- spectrum(coweeta$er.percent,log="no")  
    plot((1/ssp$freq),ssp$spec,type="l")
    points((1/ssp$freq),ssp$spec)
    abline(v=365)
    ssp$period<-1/ssp$freq
    
    ssp.df<-data.frame("period"=ssp$period,"spec"=ssp$spec)
    ssp.df[which.max(ssp.df$spec),]  # 375 days
    
    coweeta.errors<-ssp.df
    
#' MARCELL EXPERIMENTAL FOREST
#' Using the small reference catchment S2, and long-term rain gage on S2 (Station South_PCP)
#' https://www.fs.usda.gov/rds/archive/products/RDS-2020-0062/_metadata_RDS-2020-0062.html
#' use record concurrent with S2 streamflow data, March 1961-present
#' 9.6ha
#' Snow-17 Parameters: nse=0.18, r2=0.37, sd of error = 208mm SWE
#' Brute force search method within the bounds on parameters taken from He et al. 2011: j.advwatres.2010.10.002
#' Parmaters were additionally bounded based on parameter constraints and physical interpretation of parameters given in the NWS Snow-17 Calibration manual (https://www.weather.gov/media/owp/oh/hrl/docs/1_Anderson_CalbManual.pdf, last accessed April 2023)
#' Anderson, E.A., 2002. Calibration of conceptual hydrologic models for use in river forecasting. Office of Hydrologic Development, US National Weather Service, Silver Spring, MD, 372.
#' 
#' Parameter bounds and sampling scheme for brute force search: 
#' "scf"=round(runif(1,0.5,1.8),2)
#' "PXTEMP"=runif(1,-2,2)
#' "MFMAX"=round(runif(1,0.8,1.2),2)
#' "MFMIN"=round(runif(1,0.1,0.6),2)
#' "UADJ"=round(runif(1,0.03,0.19),2)
#' "MBASE"=0
#' "TIPM"=round(runif(1,0.05,0.15),2)
#' "PLWHC"=round(runif(1,0.015,0.055),2)
#' "NMF"=max(c(rnorm(1,0.15,0.05),0.05))
#' "DAYGM"=round(runif(1,0,0.3),2)
#' 
#' 56,010 parameter combinations attempted
#'       scf    PXTEMP     MFMAX     MFMIN      UADJ     MBASE      TIPM     PLWHC       NMF     DAYGM 
#'    1.7300000 1.0996366 1.0200000 0.2800000 0.0900000 0.0000000 0.0700000 0.0500000 0.2051719 0.0000000 
    
  # read in testing and training data: 
  # This includes columns for date, daily streamflow (mm, labeled Qmm), daily precipitation (mm, labeled Pmm), and mean daily temperature (degC, labeled tmean)
  # For MEF, an additional column is RAIM: Rain+Melt, this is modeled output from a brute-force calibrated Snow17 model, using the temperature forcings, and measured against the aggregated SWE on the S2 catchment. 
    load("Data/MEFtest.RData") # called mef.test
    load("Data/MEFtrain.RData") # called mef.train
  
  # create the Takens vectors: 
    m<-5  # 3 days of streamflow embedding dim
    td<-1 # 1 day time lag
    days.precip<-7  # 7 days of past precip. Will also add one day of perfectly known QPF. 
    days.temp<-5  # 3 days of past observed daily mean temperature.
    
    takSamp<-buildTakens(mef.train$Qmm, embedding.dim = m, time.lag = td)  # the training attractor for streamflow
    takPred<-buildTakens(mef.test$Qmm, embedding.dim = m, time.lag = td)  # the testing attractor for streamflow
    
    precipTakTrain<-buildTakens(mef.train$raim,embedding.dim=days.precip,time.lag=1)  # precip training vectors
    precipTakPred<-buildTakens(mef.test$raim, embedding.dim = days.precip, time.lag = 1)  # precip testing vectors
    
    # 1 day of QPF: 
    precipTakTrain<-cbind(precipTakTrain,c(precipTakTrain[2:length(precipTakTrain[,1]),days.precip],10))
    precipTakPred<-cbind(precipTakPred,c(precipTakPred[2:length(precipTakPred[,1]),days.precip],10))
    
    # now the temperature: 
    tempTrain<-buildTakens(mef.train$tmean, embedding.dim = days.temp, time.lag = td)
    tempTest<-buildTakens(mef.test$tmean, embedding.dim = days.temp, time.lag = td)
    
    # align the vectors temporally: When days.precip is longer than the others. 
    tempTrain<-tempTrain[(days.precip-(days.temp-1)):(length(tempTrain[,1])),] 
    tempTest<-tempTest[(days.precip-(days.temp-1)):(length(tempTest[,1])),] 
    takSamp<-takSamp[(days.precip-(m-1)):(length(takSamp[,1])),] 
    takPred<-takPred[(days.precip-(m-1)):(length(takPred[,1])),] 
    
  # Implement the algorithm: 
    longterm.f<-longterm.forecaster(takSamp,takPred,precipTakTrain,precipTakPred,tempTest,tempTrain)
  
  # Fit statistics:
    nse.fcn(mef.test$Qmm[8:(length(longterm.f)+2)],longterm.f[6:length(longterm.f)])  
    r2.fcn(mef.test$Qmm[8:(length(longterm.f)+2)],longterm.f[6:length(longterm.f)])
    
  # do spectral analysis: 
    error<-(longterm.f[6:length(longterm.f)]-mef.test$Qmm[8:(length(longterm.f)+2)])
    plot(error,type="l",ylim=c(-10,10))  # need to set limit because for the very large streamflow event, we WAY undersimulate. 
    
    ssp <- spectrum(error,log="no")   # do the spectrum on the raw error bc percentage is not available (many zeroes so get divide-by-zeroes)
    plot((1/ssp$freq),ssp$spec,type="l")
    points((1/ssp$freq),ssp$spec)
    abline(v=365)
    abline(v=365/2)
    abline(v=365/4)
    abline(v=365*2)
    
    ssp$period<-1/ssp$freq
    
    ssp.df<-data.frame("period"=ssp$period,"spec"=ssp$spec)
    ssp.df[which.max(ssp.df$spec),]  # 375 days
    head(ssp.df[order(ssp.df$spec,decreasing=TRUE),],5)
    
    marcell.errors<-ssp.df
    
    # Plot streamflow data: 
    marcell<-data.frame("date"=mef.test$date[8:(length(longterm.f)+2)],"obs"=mef.test$Qmm[8:(length(longterm.f)+2)],"sim"=longterm.f[6:length(longterm.f)])
    
    mef00<-subset(marcell,year(date)>=2000 & year(date)<=2001)   # only plot 2000 
    
    long.forPlot<-pivot_longer(mef00,cols=c("obs","sim"),names_to="Timeseries",values_to="Streamflow")
    
    mef00<-ggplot(data=long.forPlot,aes(date,Streamflow,color=Timeseries))+geom_line(lwd=0.75)+
      theme_bw(base_size=(size=20))+
      #ylim(0,20)+
      xlab("Date")+
      ylab("Streamflow [mm/day]")+
      scale_color_manual(values=c("black","red"))
    mef00
    
    png(filename="marcell_SF.png",width=10,height=6,res=300,units="in")
    mef00
    dev.off()
    