
# purpose is scaling diagrams from new created nc files with data from masks
# create Rdata files with necessay information for further processing
#

# some useful things


iepsout = 5
dot <- "."
dash <- "-"
ldash <- "_"


Rlibdir <- .libPaths()
library("lubridate", lib.loc=Rlibdir)
#library("ismev", lib.loc=Rlibdir)
#library("evd", lib.loc=Rlibdir)
library("PCICt", lib.loc=Rlibdir)
library("abind", lib.loc=Rlibdir)
#library("pryr", lib.loc=Rlibdir)
library("rlist",  lib.loc=Rlibdir)




gl.count <- function(xdata, xthres=xthres) {
  xr=length(xdata[xdata>xthres])
  return(xr)
}



gl.ReadFromNetcdf1D <- function(ncinfile, varname,maskfile = "no", maskname = "mask") {

  ncid <- nc_open(ncinfile,write=FALSE)
  var  <- ncvar_get(ncid,varid=varname)
  dimvar <- dim(var)
  ndim <- length(dimvar)
  dimvarNoTime <- dimvar[-ndim]

  if ( ndim==2 ) {
    print("1D file found")
  } else {
    print("dimensions of input file not correct; only time,stationid format allowed")
  }

  if (maskfile != "no") {

    print(maskfile)
    ncid2 <- nc_open(maskfile,write=FALSE)
    mask  <- ncvar_get(ncid2,varid=maskname)
    dimmask <- dim(mask)

    if (all(dimvarNoTime == dimmask )) { # check whether dimension match

      lmask <- mask == 1
      print(paste("number of gridpoints selected =",sum(mask)))
      var <- var[lmask,]

    } else {
      print("mask dimension and variable does not match")
      print(dimvarNoTime)
      print(dimmask)
    }


  } else {

   # data is already read in
    print("no maskfile found; output all data")

  }

  time <- ncvar_get(ncid,varid="time")
  ntime <- length(time)
  timeatt <- ncatt_get(ncid, "time")
  timeR <- gl.time2P(time,timeatt)
  sout <- list("data"=var,"time"=timeR)
  return(sout)

  }
  
  


gl.time2P <- function(time=time,timeatt=timeatt) {

timeunits <- unlist(strsplit(timeatt$units," "))

if ( timeatt$calendar == "365_day" | timeatt$calendar == "360_day" ) {
  print("use PCICt ")
  tstart <- as.PCICt(paste(timeunits[3],timeunits[4]),cal=timeatt$calendar,tz="UTC")
  if (timeunits[1] == "hours" ) {
    timeR <- as.PCICt(time*3600,origin=tstart,cal=timeatt$calendar, tz="UTC")
  } else if (timeunits[1] == "days" ) {
    timeR <- as.PCICt(time*3600*24,origin=tstart,cal=timeatt$calendar, tz="UTC")
  } else if (timeunits[1] == "seconds" ) {
    timeR <- as.PCICt(time,origin=tstart,cal=timeatt$calendar, tz="UTC")
  }

} else {

  print("use POSIXt tz=UTC")
  tstart <- as.POSIXct(paste(timeunits[3],timeunits[4]),tz="UTC")
  if (timeunits[1] == "hours" ) {
    timeR <- as.POSIXct(time*3600,origin=tstart,tz="UTC")
  } else if (timeunits[1] == "days" ) {
    timeR <- as.POSIXct(time*3600*24,origin=tstart,tz="UTC")
  } else if (timeunits[1] == "seconds" ) {
    timeR <- as.POSIXct(time,origin=tstart,tz="UTC")
  }
}

if ( timeatt$calendar == "365_day" ) {
    print("use PCICt; special with 365_day, do not use tz = UTC, looks like a bug ")
    tstart <- as.PCICt(paste(timeunits[3],timeunits[4]),cal=timeatt$calendar)
    if (timeunits[1] == "hours" ) {
      timeR <- as.PCICt(time*3600,origin=tstart,cal=timeatt$calendar)
    } else if (timeunits[1] == "days" ) {
      timeR <- as.PCICt(time*3600*24,origin=tstart,cal=timeatt$calendar)
    } else if (timeunits[1] == "seconds" ) {
      timeR <- as.PCICt(time,origin=tstart,cal=timeatt$calendar)
    }

}

return(timeR)


}




gl.monthlist <- function(cperiod=cperiod) {
  mm <- 1:12
  cperiod <- tolower(cperiod)
  if (cperiod == "all") {ltsel <-mm > 0}
  else if (cperiod == "jja" )   {ltsel <-  mm >=  6 & mm <=8  }
  else if (cperiod == "jjas")   {ltsel <-  mm >=  6 & mm <=9  }
  else if (cperiod == "djf" )   {ltsel <-  mm <=  2 | mm >=12 }
  else if (cperiod == "shy" )   {ltsel <-  mm >=  5 & mm <=10 }
  else if (cperiod == "why" )   {ltsel <-  mm <=  4 | mm >=11 }
  else if (cperiod == "son" )   {ltsel <-  mm <= 11 & mm >=9  }
  else if (cperiod == "mam" )   {ltsel <-  mm <=  5 & mm >=3  }
  else if (cperiod == "mamj")   {ltsel <-  mm <=  6 & mm >=3  }
  else if (cperiod == "amj")   {ltsel <-  mm <=  6 & mm >=4  }
  else if (cperiod == "jaso")   {ltsel <-  mm <= 10 & mm >=7  }
  else if (cperiod == "aso")   {ltsel <-  mm <= 10 & mm >=8  }
  else if (cperiod == "ndjf")   {ltsel <-  mm <=  2 | mm >=11 }
  else if (cperiod == "mjjaso") {ltsel <-  mm <= 10 & mm >=5  }
  else if (cperiod == "mjjason") {ltsel <-  mm <= 11 & mm >=5  }
  else if (cperiod == "mjjas")  {ltsel <-  mm <= 9 & mm >=5  }
  else if (cperiod == "no_djf")  {ltsel <-  mm >   2 & mm < 12 }
  else if (cperiod == "ondjfm")  {ltsel <-  mm >=   10 | mm < 4 }
  else if (cperiod == "amjjas")  {ltsel <-  mm >=   4 & mm < 10 }
  else if (cperiod == "nodjf")  {ltsel <-  mm >   2 & mm < 12 }
  else if (cperiod == "jas")  {ltsel <-  mm >   6 & mm < 10 }
  else if (cperiod == "mjj")  {ltsel <-  mm >   4 & mm < 8 }
  else if (cperiod == "feb")  {ltsel <-  mm == 2}
  else if (cperiod == "mar")  {ltsel <-  mm == 3}
  else if (cperiod == "apr")  {ltsel <-  mm == 4}
  else if (cperiod == "may")  {ltsel <-  mm == 5}
  else if (cperiod == "jun")  {ltsel <-  mm == 6}
  else if (cperiod == "jul")  {ltsel <-  mm == 7}
  else if (cperiod == "aug")  {ltsel <-  mm == 8}
  else if (cperiod == "sep")  {ltsel <-  mm == 9}
  else if (cperiod == "oct")  {ltsel <-  mm == 10}
  else if (cperiod == "nov")  {ltsel <-  mm == 11}
  else if (cperiod == "dec")  {ltsel <-  mm == 12}
  return(mm[ltsel])
}


gl.TimeMatchClosestN <- function(time1 ,time2,  hdiff=0, tdiff = 12,  iopt = 6, nsearch = 24, rmnonmatched = TRUE){
  
  # look for closest time to hdiff before
  # option 1: slow and safe
  # option 6, standard, (6 is used for legacy reasons)
  # hdiff, postive values mean shifted to earlier time in hour
   
  ntime  <- length(time1)
  ntime2 <- length(time2)
  iTimeTrans <- array(NA,dim=ntime)
  
  time1 <- time1 - hours(hdiff)   # target hours hdiff earlier
  
  if (iopt == 1) {
    for (i in 1:ntime){
      iTimeTrans[i] <- which.min(abs(time1[i]-time2))
    }
  } else if (iopt == 6 ) {
    
    j <- which.min(abs(time1[1]-time2))
    iTimeTrans[1] <- j
    js <- j
    i <- 2
    while (i <= ntime) {
      jsp <- min(js + nsearch,ntime2)
      jsmall <- which.min(abs(time1[i]-time2[js:jsp]))
      if (difftime(time1[i],time2[js-1+jsmall],units="hours") > tdiff) {
        jsmall <- which.min(abs(time1[i]-time2[js:ntime2]))
      }
      js <-   js - 1 + jsmall
      iTimeTrans[i] <- js    #
      i=i+1
      #      if (i %% 50000 == 0 ) print(paste("time step = ",i," from ",ntime))
    }
 }
  
  time1_test <- time2[iTimeTrans[1:ntime]]
  hourdiff <- difftime(time1_test,time1,units="hours")
  ltsel <- abs(hourdiff) > tdiff
  
  if (any(ltsel)){
    print(paste(hourdiff[ltsel],timeH_test[ltsel],timeH[ltsel]))
    if ( rmnonmatched) {print("nonmatching removed")
      iTimeTrans[ltsel] <- NA
      time1_test <- time2[iTimeTrans[1:ntime]]
      hourdiff <- difftime(time1_test,time1,units="hours")}
  }
  
  
  print(paste("mean = ",formatC(mean(hourdiff,na.rm=TRUE),format="f",digits =3), "hours after target hour" ))
  print(paste("max  = ",formatC(max(hourdiff,na.rm=TRUE),format="f",digits =3),  "hours after target hour " ))
  print(paste("min  = ",formatC(min(hourdiff,na.rm=TRUE),format="f",digits =3),  "hours after target hour" ))
  if (tdiff == 0) {print("exact time matching; are you sure?")}
  
  return(iTimeTrans)
  
}



gl.pr_extremes <- function(precip,pthres = 0.1, nsample = 1, nseed = -999, nblock = 1, indtime = 2, narm=TRUE) {

  zpctl_fac <- 10**0.1
  pctls_extremes <- as.double(1e0-(1e0/zpctl_fac)**seq(1,75,1))
  npctls_extremes <- length(pctls_extremes)
  pvar_extremes_boot <-  array(NA,c(npctls_extremes,nsample))
  pvar_extremes_boot_cond <-  array(NA,c(npctls_extremes,nsample))
  pthres_freq_boot <- array(NA,nsample)

  print("this only works if variable is 2D, with nblocks in indtime index; default = 2 !!!!!!")
  if (narm) { precip[precip < -10] <- NA }
  
# do a bootstrap resample of
  if ( nsample > 1 ){ 
    if (nseed != -999) {set.seed(nseed)}
    ntime = dim(precip)[indtime] # denk 2
    print(paste("number of time steps to be resampled =",ntime))

    for (isample in 1:nsample ){

       if (isample%%10 == 0) {
		   print(paste("sample ",isample," from ", nsample, "memory used = "))
           print(mem_used())
           }
           

       if (nblock > 1) {
         nrb <- ntime/nblock
         iper <- rep(nblock*sample(1:nrb-1,nrb,replace=TRUE),each=nblock) + rep(1:nblock)
         if (isample == 1) {print("block resampling")}
       } else {
         if (isample == 1) {print("no block resampling")}
         iper <- sample(1:ntime,ntime,replace = TRUE)
       }

# dit lijkt fout te zijn, iper in tweede dimensie

       if (indtime == 1) {pr_resample <- precip[iper,]}
       if (indtime == 2) {pr_resample <- precip[,iper]}
       prc <- pr_resample[!is.na(pr_resample)]
       nall <- length(prc)
     
       pvar_extremes <- quantile(prc, probs = pctls_extremes, type=5,names=F)
        for (ip in 1:npctls_extremes) {
         nmin = ceiling(1/(1-pctls_extremes[ip]))
         if (nall < 2*nmin) {pvar_extremes[ip] <- NA}
       }
       pvar_extremes_boot[,isample] <- pvar_extremes
       
       prc <- prc[prc>pthres]
       nwet <- length(prc)
       pthres_freq_boot[isample] <- nwet/nall
 
       pvar_extremes <- quantile(prc, probs = pctls_extremes, type=5,names=F)
        for (ip in 1:npctls_extremes) {
         nmin = ceiling(1/(1-pctls_extremes[ip]))
         if (nwet < 2*nmin) {pvar_extremes[ip] <- NA}
       }
       pvar_extremes_boot_cond[,isample] <- pvar_extremes
      

     }
     rm(pr_resample)
     rm(prc)


  }

  print("not bootstrapped")
  
  prc <- precip[!is.na(precip)]
  pvar_extremes <- quantile(prc, probs = pctls_extremes, type=5,names=F)
  for (ip in 1:npctls_extremes) {
   nmin = ceiling(1/(1-pctls_extremes[ip]))
   if (length(prc) < 2*nmin) {pvar_extremes[ip] <- NA}
  }

  nall <- length(prc)
  prc <- prc[prc>pthres]
  nwet <- length(prc)

  pvar_extremes_cond <- quantile(prc, probs = pctls_extremes, type=5,names=F)
  for (ip in 1:npctls_extremes) {
   nmin = ceiling(1/(1-pctls_extremes[ip]))
   if (length(prc) < 2*nmin) {pvar_extremes_cond[ip] <- NA}
  }

  rm(prc)
  ExtremesOut <- list("pctls_extremes" =pctls_extremes,
                      "pvar_extremes" = pvar_extremes,
                      "pvar_extremes_boot" = pvar_extremes_boot,
                      "pvar_extremes_boot_cond" = pvar_extremes_boot_cond,
                      "pthres_freq_boot" = pthres_freq_boot,
                      "pvar_extremes_cond" = pvar_extremes_cond,
                      "nall" = nall, "nwet" = nwet, "fwet" = nwet/nall ,
                      "nblock" = nblock
                      )

  return(ExtremesOut)

}



#===============================================================================
gl.PlotScaling <- function(ScalingObj, FitObj="no",
                          xlimits = c(4,24),
                          ylimits = c(0.6,120),
                          fitopt = c(TRUE,TRUE,3),
                          plotfile = "ScalingPlotDummy.eps",
                          plottitlelt = c(" "," "),
                          lsfactor = 1,
                          loessrange = c(-999,-999),
                          ylabel = "hourly precipitation [mm/hour]" ,
                          xlabel = "dew point temperature [degC]",
                          rlspan = 0.75) {


  nbins = ScalingObj$nbins
  xvar_mean = ScalingObj$xvar_mean
  pr_quant = ScalingObj$pr_quant
  NFreqExc = ScalingObj$NFreqExc
  ExcTres = ScalingObj$ExcTres
  temp_quant = ScalingObj$temp_quant
  temp_paired_range = ScalingObj$temp_paired_range
  temp_paired_wet_range = ScalingObj$temp_paired_wet_range
  temp_paired_Pgt1_range = ScalingObj$temp_paired_Pgt1_range
  temp_paired_Pgt5_range = ScalingObj$temp_paired_Pgt5_range
  temp_paired_Pgt10_range = ScalingObj$temp_paired_Pgt10_range
  temp_paired_Pgt20_range = ScalingObj$temp_paired_Pgt20_range
  temp_paired_Pgt40_range = ScalingObj$temp_paired_Pgt40_range
  pquant = ScalingObj$pquant

  fileout <- plotfile
  source("initeps.R",local=TRUE)
  options(scipen=2)
  par(cex=1.4)

  ylimits <- ylimits*lsfactor

  plot(xvar_mean,pr_quant[,1],log ="y",xlim=xlimits,ylim=ylimits,type="p",col="cyan",
       ylab = ylabel, xlab=xlabel, xaxs="i",pty="s",pch=16,yaxs="i")

  xcommon <- seq(-10,30,0.1)
  CCrate <- 0.065
  yscaling7 <- lsfactor*2*(1+CCrate)**xcommon # plot some reference lines CC and 2CC
  yscaling14 <- lsfactor*(1+2*CCrate)**xcommon
  lines(xcommon,yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,2*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,4*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,8*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,2*yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,4*yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,8*yscaling7,lty=3,lwd=0.5,col="black")

  abline(v=10,col="grey",lty=3)
  abline(v=15,col="grey",lty=3)
  abline(v=20,col="grey",lty=3)


  points(xvar_mean,pr_quant[,2],col="blue",pty="s",pch=16)
  points(xvar_mean,pr_quant[,3],col="magenta",pty="s",pch=16)

  pcha <- c(16, 18,16, 15, 16,18,16 )
  cexa <- c(0.0, 1.2,1.5, 1.8, 1.5,1.2 ,0.0)

  yline <- 1.1*rep(ylimits[1], length(temp_paired_range))
  lines(temp_paired_range,yline,col="black",lwd=2)
  points(temp_paired_range,yline,pch = pcha,
         col = "black",cex =  0.6*cexa)

  pyfac <- 1.14
  yline <- pyfac*yline
  lines(temp_paired_wet_range,yline,col="blue",lwd=2)
  points(temp_paired_wet_range,yline,pch = pcha,
         col = "blue",cex = 0.6*cexa)
         text(temp_paired_wet_range[1],yline[1],">0.1",pos=2,cex=0.7)

  yline <- pyfac*yline
  lines(temp_paired_Pgt5_range,yline,col="orange",lwd=2)
  points(temp_paired_Pgt5_range,yline,pch = pcha,
         col = "orange",cex =  0.6*cexa)
         text(temp_paired_Pgt5_range[1],yline[1],">5",pos=2,cex=0.7)

  yline <- pyfac*yline
  lines(temp_paired_Pgt20_range,yline,col="red",lwd=2)
  points(temp_paired_Pgt20_range,yline,pch = pcha,
         col = "red",cex =  0.6*cexa)

  text(temp_paired_Pgt20_range[1],yline[1],">20",pos=2,cex=0.7)

  yline <- pyfac*yline
  lines(temp_paired_Pgt40_range,yline,col="magenta",lwd=2)
  points(temp_paired_Pgt40_range,yline,pch = pcha,
         col = "magenta",cex =  0.6*cexa)
         text(temp_paired_Pgt40_range[1],yline[1],">40",pos=2,cex=0.7)

  print(paste("rlspan =",rlspan))
  if (rlspan > 0 ) {
  pcols <- c("cyan","blue","magenta")
  for (ii in 1:3) {
    ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > temp_paired_Pgt5_range[1] & xvar_mean <temp_paired_Pgt5_range[7]
    if (lsfactor < 0.5)
       {ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > temp_paired_Pgt1_range[1] & xvar_mean <temp_paired_Pgt1_range[7]}
    if (loessrange[1] != -999 & loessrange[2] - loessrange[1] > 5  )
      {ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > loessrange[1] & xvar_mean <loessrange[2]}
      
    if (any(ltsel)) { 
		if  (range(xvar_mean[ltsel])[2] - range(xvar_mean[ltsel])[1] > 5) {
			y.loess <- loess(pr_quant[ltsel,ii]~xvar_mean[ltsel], span=rlspan, data.frame(x=xvar_mean[ltsel], y=pr_quant[ltsel,ii]))
			y.predict <- predict(y.loess, data.frame(x=xvar_mean[ltsel]))
			lines(xvar_mean[ltsel],y.predict,col=pcols[ii],lwd=3) }
    }
  }
  }

  lfit = fitopt[1]
  lweightfit = fitopt[2]
  ifits = fitopt[3:length(fitopt)]

  # and plots some fits
  if (! is.list(FitObj) | lfit == FALSE ){
    print("no fitting data found or fitting set to false ")
  } else {

    print("FITTING DATA FOUND")
    Pmin_lm <- FitObj$Pmin_lm
    Pmax_lm <- FitObj$Pmax_lm
    Tmin_lm <- FitObj$Tmin_lm
    Tmax_lm <- FitObj$Tmax_lm

    if (lweightfit){
      Pmin_lm <- FitObj$Pmin_lm_weight
      Pmax_lm <- FitObj$Pmax_lm_weight
      Tmin_lm <- FitObj$Tmin_lm_weight
      Tmax_lm <- FitObj$Tmax_lm_weight
    }

    for (ifit in ifits){
    for (ii in 1:3) {
      xtmp <- c(Tmin_lm[ii,ifit],Tmax_lm[ii,ifit])
      ytmp <- c(Pmin_lm[ii,ifit],Pmax_lm[ii,ifit])
      lines(xtmp,ytmp,col = pcols[ii],cex =  0.7*cexa, lty =2)
      print(xtmp)
      print(ytmp)
    }
    }


  } #endif fitting data

  print("plotting title")
  text(xlimits[1]+0.5,0.80*ylimits[2],plottitlelt[1],pos=4,cex=1.1)
  text(xlimits[1]+0.5,0.60*ylimits[2],plottitlelt[2],pos=4,cex=1.1)

  dev.off()

}



#=======================================================================================
gl.Scaling <- function(bindef=0, pquant = c(0.9,0.99,0.999,0.9999), temp_paired  ,
                       pr_paired  , pthres = 0.1 ) {


  if (bindef[1] == 0){
    bindef <- list("bin_width" = 2, "bin_min" = -10, "bin_max" = 30, "bin_step" = 1 )
  }
  # first do bin definition
  bins_lowbound <- seq(bindef$bin_min-bindef$bin_width/2,
                        bindef$bin_max-bindef$bin_width/2,bindef$bin_step)
  bins_upbound  <- bins_lowbound + bindef$bin_width
  nbins <- length(bins_lowbound)

  ExcTres <- c(10,20,40,60)
  NFreqExc <- array(NA,c(nbins,length(ExcTres)))

  nquant <- length(pquant)
  temp_quant    <- array(NA,c(nbins,3))
  ndata_bin_all <- array(NA,nbins)
  ndata_bin_wet <- array(NA,nbins)
  pr_quant   <- array(NA,c(nbins,nquant))

  for (ib in 1:nbins) {

   # select data
   iisel <- temp_paired > bins_lowbound[ib] & temp_paired <= bins_upbound[ib]
   iiselw <- iisel & pr_paired >= pthres

   ndata_bin_all[ib] <- sum(iisel, na.rm = TRUE)
   ndata_bin_wet[ib] <- sum(iiselw, na.rm = TRUE)

   pr_sel   <-   pr_paired[iiselw]
   temp_sel <- temp_paired[iiselw]

   pr_quant[ib,]   <-   quantile(  pr_sel,probs=pquant,type=5,na.rm=TRUE)
   temp_quant[ib,] <-   quantile(temp_sel,probs=c(0.05,0.5,0.95),type=5,na.rm=TRUE)
   xvar_mean <- temp_quant[,2]

   # apparent the quantile function always outputs a percentile, set to NA if not sufficient values
   for (ip in 1:length(pquant)){
     nmin = 2*ceiling(1/(1-pquant[ip]))
     if (ndata_bin_wet[ib]<nmin) {pr_quant[ib,ip] <- NA}
   }

   for (it in 1:length(ExcTres)) {
     NFreqExc[ib,it] <- gl.count(pr_sel,ExcTres[it]) # / length(pr_sel)
   }
 }

 tquant = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
 temp_paired_range       <-quantile(temp_paired,probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_wet_range   <-quantile(temp_paired[pr_paired >= pthres],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_Pgt1_range  <-quantile(temp_paired[pr_paired > 1.],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_Pgt5_range  <-quantile(temp_paired[pr_paired > 5.],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_Pgt10_range <-quantile(temp_paired[pr_paired > 10.],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_Pgt20_range <-quantile(temp_paired[pr_paired > 20.],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_Pgt40_range <-quantile(temp_paired[pr_paired > 40.],probs = tquant,type = 5,na.rm = TRUE)

 pr_PCTL     <- quantile(pr_paired,                      probs = c(0.5,0.9,0.99,0.999), type = 5, na.rm = TRUE)
 pr_PCTL_wet <- quantile(pr_paired[pr_paired >= pthres], probs = c(0.5,0.9,0.99,0.999), type = 5, na.rm = TRUE)
 temp_paired_PCTL1_range <-quantile(temp_paired[pr_paired > pr_PCTL[1]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL2_range <-quantile(temp_paired[pr_paired > pr_PCTL[2]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL3_range <-quantile(temp_paired[pr_paired > pr_PCTL[3]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL4_range <-quantile(temp_paired[pr_paired > pr_PCTL[4]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL1_wet_range <-quantile(temp_paired[pr_paired > pr_PCTL_wet[1]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL2_wet_range <-quantile(temp_paired[pr_paired > pr_PCTL_wet[2]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL3_wet_range <-quantile(temp_paired[pr_paired > pr_PCTL_wet[3]],probs = tquant,type = 5,na.rm = TRUE)
 temp_paired_PCTL4_wet_range <-quantile(temp_paired[pr_paired > pr_PCTL_wet[4]],probs = tquant,type = 5,na.rm = TRUE)

 ScalingOut <- list("bins_lowbound" = bins_lowbound, "bins_upbound" = bins_upbound, "bins" = nbins )
 ScalingOut$xvar_mean = xvar_mean
 ScalingOut$pr_quant = pr_quant
 ScalingOut$NFreqExc = NFreqExc
 ScalingOut$ExcTres = ExcTres
 ScalingOut$temp_quant = temp_quant
 ScalingOut$temp_paired_range = temp_paired_range
 ScalingOut$temp_paired_wet_range = temp_paired_wet_range
 ScalingOut$temp_paired_Pgt1_range = temp_paired_Pgt1_range
 ScalingOut$temp_paired_Pgt5_range = temp_paired_Pgt5_range
 ScalingOut$temp_paired_Pgt10_range = temp_paired_Pgt10_range
 ScalingOut$temp_paired_Pgt20_range = temp_paired_Pgt20_range
 ScalingOut$temp_paired_Pgt40_range = temp_paired_Pgt40_range
 ScalingOut$pquant = pquant
 ScalingOut$PthresFreq <- ndata_bin_wet/ndata_bin_all
 ScalingOut$Nall <- ndata_bin_all
 ScalingOut$Nwet <- ndata_bin_wet

 ScalingOut$pr_PCTL <- pr_PCTL
 ScalingOut$pr_PCTL_wet <- pr_PCTL_wet
 ScalingOut$temp_paired_PCTL1_range  <- temp_paired_PCTL1_range
 ScalingOut$temp_paired_PCTL2_range  <- temp_paired_PCTL2_range
 ScalingOut$temp_paired_PCTL3_range  <- temp_paired_PCTL3_range
 ScalingOut$temp_paired_PCTL4_range  <- temp_paired_PCTL4_range
 ScalingOut$temp_paired_PCTL1_wet_range <- temp_paired_PCTL1_wet_range
 ScalingOut$temp_paired_PCTL2_wet_range <- temp_paired_PCTL2_wet_range
 ScalingOut$temp_paired_PCTL3_wet_range <- temp_paired_PCTL3_wet_range
 ScalingOut$temp_paired_PCTL4_wet_range <- temp_paired_PCTL4_wet_range

# returns median of temperature range of all events exceeding different P threshold based
# on absolute percentiles, same information as above but condensed

 prcTthres <- c(NA,pthres,pr_PCTL)
 names(prcTthres) <- c("all","wet",names(pr_PCTL))
 TdStat    <- array(NA,length(prcTthres))
 TdStat[1] <- temp_paired_range[4]
 TdStat[2] <- temp_paired_wet_range[4]
 TdStat[3] <- temp_paired_PCTL1_range[4]
 TdStat[4] <- temp_paired_PCTL2_range[4]
 TdStat[5] <- temp_paired_PCTL3_range[4]
 TdStat[6] <- temp_paired_PCTL4_range[4]

 ScalingOut$prcThres <- prcTthres
 ScalingOut$TdStat   <- TdStat

 return(ScalingOut)

}



#=======================================================================================
gl.PDF <- function      (varin ,
                         ntoint = 10,
                         ymax = 100,
                         ymin = 0,
                         xcoarse_in = c(seq(0.0,1,0.1),1+seq(1,40,1)**2*0.1),
                         undef = -999
                     ) {

  # check things

  print("HistFine counts between >=xFine[i] & <xFine[i+1]")


  # remove NA and undefine from data from dataid
  print(paste("removed NA and undef = ",undef,"from data"))
  nrint = 100000
  varin <- round(varin*nrint)/nrint
  varin <- varin[! is.na(varin) & varin != undef]

  if (min(varin,na.rm=TRUE) < ymin ) {
    print(paste("WARNING: minimum data =",min(varin,na.rm=TRUE)," and ymin =", ymin))}

  if (max(varin,na.rm=TRUE) > ymax ) {
    print(paste("WARNING: maximum data =",max(varin,na.rm=TRUE)," and ymax =", ymax))}



  # notoint gives steps for fine resolution bins, 100 means steps of 0.01
  nymin <- floor(ymin*ntoint)
  nymax <- floor(ymax*ntoint)
  xfine <- round(seq(nymin,nymax,1))
  varint <- as.vector(floor(varin*ntoint))

  # built histogram for xfine, remove data outside ymin and ymax
  lsel <- varint < nymin | varint > nymax
  if (any(lsel))
       {print(paste("number of data outside boundaries =",length(which(lsel))))
        print("DO YOU WANT THIS, otherwise updata ymin and/or ymax ?")}

  xtable <- rle(sort(varint[!lsel]))
  histtable <- array(0,length(xfine))
  histtable[match(xtable$values,xfine)]<- xtable$lengths
  nFine <- sum(histtable)

 # do coarse graining of histogram

  xx <- xcoarse_in[xcoarse_in >= ymin & xcoarse_in <= ymax ]
  nbins <- length(xx)
  xx <- c(-9999,xx) # extend by one
  xx <- floor(ntoint*xx+1e-3)  # add a little bit to avoid rounding errors
  print(xx)

  HistCoarse <- array(0,nbins)
  xcoarse <- array(NA,nbins)
  ywid  <- array(NA,nbins)

  for (ip in 1:nbins){
    ntest <- which(xfine > xx[ip] & xfine <= xx[ip+1] )
    if (length(ntest)>0){
      HistCoarse[ip] <- sum(histtable[ntest])
      xcoarse[ip]    <- median(xfine[ntest])
      ywid[ip]       <- length(unique(xfine[ntest])) }
    }

  print(xcoarse)

  nCoarse <- sum(HistCoarse)
  PDFcoarse <- ntoint*HistCoarse/(ywid*nCoarse)

  xcoarse <- xcoarse / (1.*ntoint) # transfer back to original values
  xfine <- xfine / (1.*ntoint)
  
  print(xcoarse)

  PDFObj <-  list("xCoarse"= xcoarse)
  PDFObj$pdfCoarse = PDFcoarse
  PDFObj$HistCoarse = HistCoarse
  PDFObj$xFine = xfine
  PDFObj$HistFine  = histtable
  PDFObj$pdfFine  = ntoint*histtable/nFine
  PDFObj$nFine <- nFine
  PDFObj$nCoarse <- nCoarse
  PDFObj$wCoarse <- ywid
  PDFObj$xtable <- xtable
  return(PDFObj)
}



gl.getAreaTitle <- function(areaname)

{
areatitle <- areaname

areas = c("BENELUXbox.land.oroglt400",      #1
          "GERLbox.land.oroglt400",         #2
          "FRANCEWESTbox.land.oroglt400",   #3
          "FRANCESOUTHbox.land.oroglt400",  #4
          "MEDSEAbox.land.oroglt400",       #5
          "HBOUND_SKIP100.land.oroglt400",  #6
          "HBOUND_SKIP100.land.oroggt400",  #7
          "MEDSEAbox.sea.all" ,             #8
          "NWEURbox.land.oroglt400" ,       #9
          "BENELUXbox.land.all"             #10        
          )

areastitles = c("BENELUX land <400m",   #1
               "GERL land <400m" , #2
               "FRANCEWEST land <400m", #3
               "FRANCESOUTH land <400m", #4
               "MEDSEA land <400m",  #5
               "ALL DOMAIN land <400m", #6
               "ALL DOMAIN land >400m" , #7
               "MEDSEA sea ", #8
               "NWEUR land <400m",   #9
               "BENELUX land" 
               )
 
 if (areaname %in% areas) {             
 
  areatitle = areastitles[areas==areaname]       
         
}
return(areatitle)

}



gl.CDFexceedancePr <- function(td,tas,pr,pctl=list(0.1,c(0.5,0.9,0.95, 0.99,0.995,0.999)),narm=FALSE)
  
{
  
  dim(pr)  <- length(pr)
  dim(tas) <- length(tas)
  dim(td)  <- length(td)
  
  if (narm) {
	   lsel <- is.na(pr) | is.na(td) | is.na(pr) | pr < -1 | td < -90 | tas < -90 
	   print(paste("removed ",sum(lsel), "from ", length(lsel),"data"))
	   pr <-    pr[!lsel]
	   tas <-  tas[!lsel]
	   td <-    td[!lsel]	   
   }   
  
  prthres <- c(-1,pctl[[1]],quantile(pr[pr>pctl[[1]]],pctl[[2]]))
  names(prthres) <- c("all","wet",paste0(pctl[[2]]*100,"%"))
  tddepres <- tas - td 
  xseq <- seq(0.001,1,0.002)
  nprthres <- length(prthres)
  tdrange <- array(NA,c(nprthres,length(xseq)))
  tasrange <- array(NA,c(nprthres,length(xseq)))
  tddepresrange <- array(NA,c(nprthres,length(xseq)))
  
  for (ii in 1:nprthres){
    
    lsel <- pr>prthres[ii]
    if (sum(lsel,na.rm=TRUE) > 100 ){
    tdrange[ii,]   <- quantile(td[lsel],xseq)
    tasrange[ii,] <- quantile(tas[lsel],xseq)
    tddepresrange[ii,]  <- quantile(tddepres[lsel],xseq)
    }
    
  }
  
  ExcStat <- list("xseq" = xseq, "tdDist" = tdrange )
  ExcStat$tasDist <- tasrange
  ExcStat$tddepresDist <-  tddepresrange
  ExcStat$prthres <- prthres
  return(ExcStat)
  
}



gl.PlotCDFexceedancePr <- function(xdata,tempdata,prthres,
                                   xtitle = "temperature [degC]",
                                   ytitle = "cumulative fraction [0:1]",
                                   xrange = NA,
                                   pnamesplot = c("all","wet", "90%","95%", "99%" ),
                                   fileout = "dummy.eps",
                                   ptitle = c("","","")
)
 
  
{  
source("initeps.R",local=TRUE)
options(scipen=2)
par(cex=1.4)

cols <- c("black",   "green",   "blue"  ,  "magenta", "red")


if (is.na(xrange[1])) {
  tmp <- tempdata[1,xdata>0.01 & xdata < 0.99]
  xlimits <- c(min(tmp,na.rm=TRUE),max(tmp,na.rm=TRUE))  
} else {
  xlimits = xrange
}
  
ii <- 1
lsel <- names(prthres) == pnamesplot[ii]
plot(tempdata[lsel,],xdata,type="l",lwd =2, col=cols[ii],
     xlab=xtitle,ylab=ytitle, xlim = xlimits)
     
      
     
plegend <- "-"           
for (ii in 2:length(pnamesplot)){
  lsel <- names(prthres) == pnamesplot[ii]
  plegend <- c(plegend,format(unname(prthres[lsel]),digits=3))
  lines(tempdata[lsel,],xdata,lwd =3, col=cols[ii])
  print(paste(ii,names(prthres[lsel]),cols[ii]))
  
}

plegend <- paste0(pnamesplot," (", plegend,")")

legend(0.99*xlimits[2],0.01, legend=plegend,
       col=cols[1:length(pnamesplot)], lty=1, cex=0.8,xjust=1,yjust=0,lwd=2,bg="white")


dx = xlimits[2]-xlimits[1]
text(xlimits[1]+0.01*dx, 0.99,ptitle[1],pos=4,cex=0.8)
text(xlimits[1]+0.01*dx, 0.95,ptitle[2],pos=4,cex=0.8)
text(xlimits[1]+0.01*dx, 0.91,ptitle[3],pos=4,cex=0.8)
grid(col="grey20")
 

dev.off()

}


gl.Plot_POE_cluster <- function( ExtremesList, lcond = TRUE, fileout = "POE_cluster.eps", 
                            cols = c( "cyan", "orange", "blue", "magenta"),
                            xrange=c(1,5e-6), yrange=c(0,80),
                            xlabel="pooled fraction of exceedance [0..1]",ylabel="hourly rainfall [mm/hour]",
                            luncertainty = FALSE)
  {
  	
  
  setEPS()
  { 
  source("initeps.R",local=TRUE)
   
     # cairo_ps(file=fileout,width=8,height=8,family="Arial")
     #  eps(file=fileout, width=8, height=8, horizontal=TRUE)
    
    if (lcond == TRUE) {
      pvar <- "pvar_extremes_cond"
    } else {
      pvar <- "pvar_extremes"
      } 
    
    options(scipen=2)
    par(cex=1.4)
    ytop = 70
    attach(ExtremesList,2)
    ExtremesObj <- get("all",2)
    detach(ExtremesList)
    attach(ExtremesObj,2)
    pctls_extremes <- get("pctls_extremes",2)
    pvar_extremes <- get(pvar,2)
    detach(ExtremesObj)
    
    xx <- 1-pctls_extremes
    plot(xx,pvar_extremes,type="l",log="x",col="grey60",lwd=2, xlim=xrange,ylim=yrange,
         xlab=xlabel,ylab=ylabel)
    grid()
    
    if (luncertainty){
	      if ("pvar_extremes_boot" %in% names(ExtremesObj)){ # test onexistence
              pvar_extremes_boot_range <- apply(ExtremesObj$pvar_extremes_boot,1,quantile,probs=c(0.1,0.5,0.9),na.rm=TRUE) 
              lines(xx,pvar_extremes_boot_range[1,],lty=2,col="grey",lwd=1)
              lines(xx,pvar_extremes_boot_range[3,],lty=2,col="grey",lwd=1)
              }
          } else {
               lines(xx,1.25*pvar_extremes,col="grey",lty=2)
               lines(xx,0.75*pvar_extremes,col="grey",lty=2) 
               }
               
    cfwet = format(100*ExtremesObj$fwet,digits=3)
    text(1,yrange[2],paste0("ALL", " whf = ",cfwet),col="grey60",pos=4)
    
    for ( ii in 1:4){
      
      clusterid <- paste0("cluster",ii) 
      attach(ExtremesList,2)
      ExtremesObj <- get(clusterid,2)
      detach(ExtremesList)
      attach(ExtremesObj,2)
      pctls_extremes <- get("pctls_extremes",2)
      pvar_extremes <- get(pvar,2)
      detach(ExtremesObj)
     
      xx <- 1-pctls_extremes
      lines(xx,pvar_extremes,col = cols[ii],lwd=2)
      cfwet = format(100*ExtremesObj$fwet,digits=3)
      dy = yrange[2]/18
      text(1,yrange[2]-(ii)*dy,paste0("cl.",ii," whf = ",cfwet),col=cols[ii],pos=4)
      
      if (luncertainty){
	      if ("pvar_extremes_boot" %in% names(ExtremesObj)){ # test onexistence
              pvar_extremes_boot_range <- apply(ExtremesObj$pvar_extremes_boot,1,quantile,probs=c(0.1,0.5,0.9),na.rm=TRUE) 
              lines(xx,pvar_extremes_boot_range[1,],lty=2,col=cols[ii],lwd=1)
              lines(xx,pvar_extremes_boot_range[3,],lty=2,col=cols[ii],lwd=1)
              }
          }

      
    }
    
    
    
    
    dev.off()}

}




#===============================================================================
gl.PlotScalingBoot <- function(ScalingBoot, FitObj="no",
                           xlimits = c(4,24),
                           ylimits = c(0.6,120),
                           fitopt = c(TRUE,TRUE,3),
                           plotfile = "ScalingPlotDummy.eps",
                           plottitlelt = c(" "," "),
                           lsfactor = 1,
                           loessrange = c(-999,-999),
                           ylabel = "hourly precipitation [mm/hour]" ,
                           xlabel = "dew point temperature [degC]",
                           rlspan = 0.75) {
  
  
#  ScalingBoot <- ScalingList # for testing
  ScalingObj = ScalingBoot$raw
  
  nbins = ScalingObj$nbins
  xvar_mean = ScalingObj$xvar_mean
  pr_quant = ScalingObj$pr_quant
  NFreqExc = ScalingObj$NFreqExc
  ExcTres = ScalingObj$ExcTres
  temp_quant = ScalingObj$temp_quant
  temp_paired_range = ScalingObj$temp_paired_range
  temp_paired_wet_range = ScalingObj$temp_paired_wet_range
  temp_paired_Pgt1_range = ScalingObj$temp_paired_Pgt1_range
  temp_paired_Pgt5_range = ScalingObj$temp_paired_Pgt5_range
  temp_paired_Pgt10_range = ScalingObj$temp_paired_Pgt10_range
  temp_paired_Pgt20_range = ScalingObj$temp_paired_Pgt20_range
  temp_paired_Pgt40_range = ScalingObj$temp_paired_Pgt40_range
  pquant = ScalingObj$pquant
  
  
  # do the uncertainty calculations from the bootstrap
  
 
  nboot <- length(ScalingBoot)-1  # first one is the raw unbootstrapped outcome
  
  for ( iboot in 1:nboot){
    
    ncb <- paste0("iboot",iboot)
    xtmp <- ScalingBoot[[iboot]]$xvar_mean
    ytmp <- ScalingBoot[[iboot]]$pr_quant
    
    if (iboot == 1){
      nx <- dim(ytmp)[1]
      nq <- dim(ytmp)[2]
      pr_quant_boot <- array(NA,dim= c(nx,nq,nboot)) 
      xvar_boot     <- array(NA,dim= c(nx,nboot)) 
    }
    pr_quant_boot[,,iboot] <- ytmp
    xvar_boot[,iboot] <- xtmp
    
  }
  
  pr_uncertainty <- apply(pr_quant_boot,FUN=quantile,probs=c(0.05,0.5,0.95),na.rm=TRUE,MARGIN = c(1,2))
  ltmp <- ! apply(! is.na(pr_quant_boot),FUN=all,MARGIN=c(1,2))
  
  for (ip in 1:3) {
    ytmp <-  pr_uncertainty[ip,,]
    ytmp[ltmp] <- NA
    pr_uncertainty[ip,,] <- ytmp
  }
  
  
  fileout <- plotfile
  source("initeps.R",local=TRUE)
  options(scipen=2)
  par(cex=1.4)
  
  ylimits <- ylimits*lsfactor
  
  plot(xvar_mean,pr_quant[,1],log ="y",xlim=xlimits,ylim=ylimits,type="p",col="cyan",
       ylab = ylabel, xlab=xlabel, xaxs="i",pty="s",pch=16,yaxs="i")
  
  xcommon <- seq(-10,30,0.1)
  CCrate <- 0.065
  yscaling7 <- lsfactor*2*(1+CCrate)**xcommon # plot some reference lines CC and 2CC
  yscaling14 <- lsfactor*(1+2*CCrate)**xcommon
  lines(xcommon,yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,2*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,4*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,8*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,2*yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,4*yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,8*yscaling7,lty=3,lwd=0.5,col="black")
  
  abline(v=10,col="grey",lty=3)
  abline(v=15,col="grey",lty=3)
  abline(v=20,col="grey",lty=3)
  
  pcha <- c(16, 18,16, 15, 16,18,16 )
  cexa <- c(0.0, 1.2,1.5, 1.8, 1.5,1.2 ,0.0)
  
  yline <- 1.1*rep(ylimits[1], length(temp_paired_range))
  lines(temp_paired_range,yline,col="black",lwd=2)
  points(temp_paired_range,yline,pch = pcha,
         col = "black",cex =  0.6*cexa)
  
  pyfac <- 1.14
  yline <- pyfac*yline
  lines(temp_paired_wet_range,yline,col="blue",lwd=2)
  points(temp_paired_wet_range,yline,pch = pcha,
         col = "blue",cex = 0.6*cexa)
  text(temp_paired_wet_range[1],yline[1],">0.1",pos=2,cex=0.7)
  
  yline <- pyfac*yline
  lines(temp_paired_Pgt5_range,yline,col="orange",lwd=2)
  points(temp_paired_Pgt5_range,yline,pch = pcha,
         col = "orange",cex =  0.6*cexa)
  text(temp_paired_Pgt5_range[1],yline[1],">5",pos=2,cex=0.7)
  
  yline <- pyfac*yline
  lines(temp_paired_Pgt20_range,yline,col="red",lwd=2)
  points(temp_paired_Pgt20_range,yline,pch = pcha,
         col = "red",cex =  0.6*cexa)
  
  text(temp_paired_Pgt20_range[1],yline[1],">20",pos=2,cex=0.7)
  
  yline <- pyfac*yline
  lines(temp_paired_Pgt40_range,yline,col="magenta",lwd=2)
  points(temp_paired_Pgt40_range,yline,pch = pcha,
         col = "magenta",cex =  0.6*cexa)
  text(temp_paired_Pgt40_range[1],yline[1],">40",pos=2,cex=0.7)
  
  print(paste("rlspan =",rlspan))
  if (rlspan > 0 ) {
    pcols <- c("cyan","blue","magenta")
    for (ii in 1:3) {
      ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > temp_paired_Pgt5_range[1] & xvar_mean <temp_paired_Pgt5_range[7]
      if (lsfactor < 0.5)
      {ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > temp_paired_Pgt1_range[1] & xvar_mean <temp_paired_Pgt1_range[7]}
      if (loessrange[1] != -999 & loessrange[2] - loessrange[1] > 5  )
      {ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > loessrange[1] & xvar_mean <loessrange[2]}
      
      if (any(ltsel)) { 
        if  (range(xvar_mean[ltsel])[2] - range(xvar_mean[ltsel])[1] > 5) {
          
          xx <-   c(xvar_mean,rev(xvar_mean))
          yy <- c(pr_uncertainty[1,,ii],rev(pr_uncertainty[3,,ii]))
        #  xx <-   c(xvar_mean[ltsel],rev(xvar_mean[ltsel]))
        #  yy <- c(pr_uncertainty[1,ltsel,ii],rev(pr_uncertainty[3,ltsel,ii]))
          
          xx <- xx[!is.na(yy)   ]
          yy <- yy[!is.na(yy)   ]
          polygon(xx,yy,col="grey90",border=NA)
          
          y.loess <- loess(pr_quant[ltsel,ii]~xvar_mean[ltsel], span=rlspan, data.frame(x=xvar_mean[ltsel], y=pr_quant[ltsel,ii]))
          y.predict <- predict(y.loess, data.frame(x=xvar_mean[ltsel]))
          lines(xvar_mean[ltsel],y.predict,col=pcols[ii],lwd=3) }
        
      }
    }
  }

  points(xvar_mean,pr_quant[,1],col="cyan",pty="s",pch=16)
  points(xvar_mean,pr_quant[,2],col="blue",pty="s",pch=16)
  points(xvar_mean,pr_quant[,3],col="magenta",pty="s",pch=16)
  
  

  
  lfit = fitopt[1]
  lweightfit = fitopt[2]
  ifits = fitopt[3:length(fitopt)]
  
  # and plots some fits
  if (! is.list(FitObj) | lfit == FALSE ){
    print("no fitting data found or fitting set to false ")
  } else {
    
    print("FITTING DATA FOUND")
    Pmin_lm <- FitObj$Pmin_lm
    Pmax_lm <- FitObj$Pmax_lm
    Tmin_lm <- FitObj$Tmin_lm
    Tmax_lm <- FitObj$Tmax_lm
    
    if (lweightfit){
      Pmin_lm <- FitObj$Pmin_lm_weight
      Pmax_lm <- FitObj$Pmax_lm_weight
      Tmin_lm <- FitObj$Tmin_lm_weight
      Tmax_lm <- FitObj$Tmax_lm_weight
    }
    
    for (ifit in ifits){
      for (ii in 1:3) {
        xtmp <- c(Tmin_lm[ii,ifit],Tmax_lm[ii,ifit])
        ytmp <- c(Pmin_lm[ii,ifit],Pmax_lm[ii,ifit])
        lines(xtmp,ytmp,col = pcols[ii],cex =  0.7*cexa, lty =2)
        print(xtmp)
        print(ytmp)
      }
    }
    
    
  } #endif fitting data
  
  print("plotting title")
  text(xlimits[1]+0.5,0.80*ylimits[2],plottitlelt[1],pos=4,cex=1.1)
  text(xlimits[1]+0.5,0.60*ylimits[2],plottitlelt[2],pos=4,cex=1.1)
  
  dev.off()
  
}




#=======================================================================================
gl.smooth <- function   (xin ,
                         yin ,
                         xout ,
                         drange = c(min(xout),max(xout)),
                         rlspan = 0.6,
                         undef = -999,
                         zweights = 1.
) {


# print(paste("rlspan =", rlspan))
xx <- xin
yy <- yin
if ( all(zweights == 1.) ) {
      y.loess <- loess(yy~xx, span=rlspan)  
  } else {
	  y.loess <- loess(yy~xx, span=rlspan, weights = zweights ) 
  }
  
pnew <- data.frame(xx = xout)
ysmooth  <- predict(y.loess, pnew)
ysmooth[xout < drange[1] | xout > drange[2]] <- NA  # puts everything outside range to NA

yyout <- list("y"=ysmooth,"x"=xout)
return(yyout)

}



#===============================================================================
gl.PlotScalingDuo <- function(ScalingObj, ScalingObj2 = "NA",
                           xlimits = c(4,24),
                           ylimits = c(0.6,120),
                           plotfile = "ScalingPlotDummyDuo.eps",
                           plottitlelt = c(" "," "),
                           plabel = "a",
                           ylabel = "hourly precipitation [mm/hour]" ,
                           xlabel = "dew point temperature [degC]" ) {



  nbins = ScalingObj$nbins
  xvar_mean = ScalingObj$xvar_mean
  pr_quant = ScalingObj$pr_quant
  NFreqExc = ScalingObj$NFreqExc
  ExcTres = ScalingObj$ExcTres
  temp_quant = ScalingObj$temp_quant
  temp_paired_range = ScalingObj$temp_paired_range
  temp_paired_wet_range = ScalingObj$temp_paired_wet_range
  temp_paired_Pgt1_range = ScalingObj$temp_paired_Pgt1_range
  temp_paired_Pgt5_range = ScalingObj$temp_paired_Pgt5_range
  temp_paired_Pgt10_range = ScalingObj$temp_paired_Pgt10_range
  temp_paired_Pgt20_range = ScalingObj$temp_paired_Pgt20_range
  temp_paired_Pgt40_range = ScalingObj$temp_paired_Pgt40_range
  pquant = ScalingObj$pquant


  #dev.off()
  fileout <- plotfile
  source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
  options(scipen=2)
  par(cex=1.4)

  pctlcols <- c("cyan","blue","magenta","red")
  pctlcols2 <- c("cyan4","blue4","magenta4","red4")

  plot(xvar_mean,pr_quant[,1],log ="y",xlim=xlimits,ylim=ylimits,type="p",col="cyan",
       ylab = ylabel, xlab = xlabel, xaxs="i",pty="s",pch=16,yaxs="i")


  CCrate <- 0.065
  xcommon <- seq(0,30,0.1)

  yscaling7 <- 2*(1+CCrate)**xcommon # plot some reference lines CC and 2CC
  yscaling14 <- (1+2*CCrate)**xcommon
  lines(xcommon,yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,2*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,4*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,8*yscaling14,lty=3,lwd=0.5,col="red")
  lines(xcommon,yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,2*yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,4*yscaling7,lty=3,lwd=0.5,col="black")
  lines(xcommon,8*yscaling7,lty=3,lwd=0.5,col="black")

  abline(v=10,col="grey",lty=3)
  abline(v=15,col="grey",lty=3)
  abline(v=20,col="grey",lty=3)


  points(xvar_mean,pr_quant[,2],col="blue",pty="s",pch=16)
  points(xvar_mean,pr_quant[,3],col="magenta",pty="s",pch=16)

  pcha <- c(16, 18,16, 15, 16,18,16 )
  cexa <- c(0.0, 1.2,1.5, 1.8, 1.5,1.2 ,0.0)

  yline <- 1.1*rep(ylimits[1], length(temp_paired_range))
  pyfac <- 1.25
  lines(temp_paired_wet_range,yline,col="blue",lwd=2)
  points(temp_paired_wet_range,yline,pch = pcha,
         col = "blue",cex = 0.6*cexa)
  text(temp_paired_wet_range[7],yline[7],">0.1",pos=4,cex=0.7)

  yline <- pyfac*yline
  lines(temp_paired_Pgt5_range,yline,col="orange",lwd=2)
  points(temp_paired_Pgt5_range,yline,pch = pcha,
         col = "orange",cex =  0.6*cexa)
  text(temp_paired_Pgt5_range[1],yline[1],">5",pos=2,cex=0.7)

  yline <- pyfac*yline
  lines(temp_paired_Pgt20_range,yline,col="red",lwd=2)
  points(temp_paired_Pgt20_range,yline,pch = pcha,
         col = "red",cex =  0.6*cexa)
  text(temp_paired_Pgt20_range[1],yline[1],">20",pos=2,cex=0.7)


  rlspan = 0.75
  pcols <- c("cyan","blue","magenta")
  for (ii in 1:3) {
    ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > temp_paired_Pgt5_range[1] & xvar_mean <temp_paired_Pgt5_range[7]
    if (sum(ltsel) >= 7 ) {
      y.loess <- loess(pr_quant[ltsel,ii]~xvar_mean[ltsel], span=rlspan, data.frame(x=xvar_mean[ltsel], y=pr_quant[ltsel,ii]))
      y.predict <- predict(y.loess, data.frame(x=xvar_mean[ltsel]))
      lines(xvar_mean[ltsel],y.predict,col=pcols[ii],lwd=3)
    }
  }

# second objectfile

  if (ScalingObj2[[1]][1] != "NA") {

  nbins = ScalingObj2$nbins
  xvar_mean = ScalingObj2$xvar_mean
  pr_quant = ScalingObj2$pr_quant
  NFreqExc = ScalingObj2$NFreqExc
  ExcTres = ScalingObj2$ExcTres
  temp_quant = ScalingObj2$temp_quant
  temp_paired_range = ScalingObj2$temp_paired_range
  temp_paired_wet_range = ScalingObj2$temp_paired_wet_range
  temp_paired_Pgt1_range = ScalingObj2$temp_paired_Pgt1_range
  temp_paired_Pgt5_range = ScalingObj2$temp_paired_Pgt5_range
  temp_paired_Pgt10_range = ScalingObj2$temp_paired_Pgt10_range
  temp_paired_Pgt20_range = ScalingObj2$temp_paired_Pgt20_range
  temp_paired_Pgt40_range = ScalingObj2$temp_paired_Pgt40_range
  pquant = ScalingObj2$pquant


  points(xvar_mean,pr_quant[,1],col=pctlcols2[1],pty="s",pch=18)
  points(xvar_mean,pr_quant[,2],col=pctlcols2[2],pty="s",pch=18)
  points(xvar_mean,pr_quant[,3],col=pctlcols2[3],pty="s",pch=18)

  pcha <- c(16, 18,16, 15, 16,18,16 )
  cexa <- c(0.0, 1.2,1.5, 1.8, 1.5,1.2 ,0.0)

  rlspan = 0.75
  pcols <- pctlcols2
  for (ii in 1:3) {
    ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > temp_paired_Pgt5_range[1] & xvar_mean <temp_paired_Pgt5_range[7]
    if (sum(ltsel) >= 7 ) {
      y.loess <- loess(pr_quant[ltsel,ii]~xvar_mean[ltsel], span=rlspan, data.frame(x=xvar_mean[ltsel], y=pr_quant[ltsel,ii]))
      y.predict <- predict(y.loess, data.frame(x=xvar_mean[ltsel]))
      lines(xvar_mean[ltsel],y.predict,col=pcols[ii],lwd=3,lty=6)
   }
  }

  yline <- 1.2*rep(ylimits[1], length(temp_paired_range))
  lines(temp_paired_wet_range,yline,col="blue4",lwd=2,lty=6)
  points(temp_paired_wet_range,yline,pch = pcha,
         col = "blue4",cex = 0.6*cexa)

  yline <- pyfac*yline
  lines(temp_paired_Pgt5_range,yline,col="orange4",lwd=2,lty=6)
  points(temp_paired_Pgt5_range,yline,pch = pcha,
         col = "orange4",cex =  0.6*cexa)

  yline <- pyfac*yline
  lines(temp_paired_Pgt20_range,yline,col="red4",lwd=2,lty=6)
  points(temp_paired_Pgt20_range,yline,pch = pcha,
         col = "red4",cex =  0.6*cexa)

  } # end second object printing


#  text(xlimits[1],0.83*ylimits[2],plottitlelt[1],pos=4,cex=1.1)
#  text(xlimits[1],0.62*ylimits[2],plottitlelt[2],pos=4,cex=1.1)

  xfac = 0.08*(xlimits[2]-xlimits[1])
  text(xlimits[1]+0.6*xfac,0.83*ylimits[2],plottitlelt[1],pos=4,cex=1.1)
  text(xlimits[1]+0.6*xfac,0.62*ylimits[2],plottitlelt[2],pos=4,cex=1.1)
  xp <- c(xlimits[1]+0.15*xfac,xlimits[1]+0.4*xfac, xlimits[1]+0.65*xfac)
  pctlcols <- c("cyan","blue","magenta","red")
  pctlcols2 <- c("cyan4","blue4","magenta4","red4")
  points(xp,rep(0.83*ylimits[2],3),pch=16,col=pctlcols)
  points(xp,rep(0.62*ylimits[2],3),pch=18,col=pctlcols2)
  text(xlimits[1]-1.5,0.83*ylimits[2],paste0(plabel,")"),pos=2,cex=1.1,xpd=NA)

  dev.off()

}



#===============================================================================
gl.PlotScalingDuoBoot <- function(ScalingBoot1, ScalingBoot2 = "NA",
                                  xlimits = c(4,24),
                                  ylimits = c(0.6,120),
                                  plotfile = "ScalingPlotDummyDuo.eps",
                                  plottitlelt = c(" "," "),
                                  plabel = "a",
                                  ylabel = "hourly precipitation [mm/hour]" ,
                                  xlabel = "dew point temperature [degC]" ,
                                  loessrange = c(-999,-999),
                                  lpoutsideonly = TRUE ) {


cScalingIn  <- c("ScalingBoot1","ScalingBoot2")
iepsout = 5
fileout <- plotfile
source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
options(scipen=2)
par(cex=1.1)


# colors of percentile and colors of lines
pctcols <- c("cyan","blue","magenta","red", "cyan4","blue4","magenta4","red4")
dim(pctcols) <- c(4,2)
lcols <- c("blue","orange","red", "blue4","orange4","red4" )
dim(lcols) <- c(3,2)
ucols <- c("blue","brown")
ucoltrans <- c(90,90)
# create empty plot
plot(0:20,0:20+NA,log ="y",xlim=xlimits,ylim=ylimits,type="p",col="cyan",
     ylab = ylabel, xlab = xlabel, xaxs="i",pty="s",pch=16,yaxs="i")


CCrate <- 0.065
xcommon <- seq(0,30,0.1)

yscaling7 <- 2*(1+CCrate)**xcommon # plot some reference lines CC and 2CC
yscaling14 <- (1+2*CCrate)**xcommon
lines(xcommon,yscaling14,lty=3,lwd=0.5,col="red")
lines(xcommon,2*yscaling14,lty=3,lwd=0.5,col="red")
lines(xcommon,4*yscaling14,lty=3,lwd=0.5,col="red")
lines(xcommon,8*yscaling14,lty=3,lwd=0.5,col="red")
lines(xcommon,yscaling7,lty=3,lwd=0.5,col="black")
lines(xcommon,2*yscaling7,lty=3,lwd=0.5,col="black")
lines(xcommon,4*yscaling7,lty=3,lwd=0.5,col="black")
lines(xcommon,8*yscaling7,lty=3,lwd=0.5,col="black")

abline(v=10,col="grey",lty=3)
abline(v=15,col="grey",lty=3)
abline(v=20,col="grey",lty=3)


for (is in 1:2){ # loop over two bootstrap Scaling objects

    ScalingBoot <-  eval(parse(text=cScalingIn[is]))

    ScalingObj <- ScalingBoot$raw
    nbins = ScalingObj$nbins
    xvar_mean = ScalingObj$xvar_mean
    pr_quant = ScalingObj$pr_quant
    NFreqExc = ScalingObj$NFreqExc
    ExcTres = ScalingObj$ExcTres
    temp_quant = ScalingObj$temp_quant
    temp_paired_range = ScalingObj$temp_paired_range
    temp_paired_wet_range = ScalingObj$temp_paired_wet_range
    temp_paired_Pgt1_range = ScalingObj$temp_paired_Pgt1_range
    temp_paired_Pgt5_range = ScalingObj$temp_paired_Pgt5_range
    temp_paired_Pgt10_range = ScalingObj$temp_paired_Pgt10_range
    temp_paired_Pgt20_range = ScalingObj$temp_paired_Pgt20_range
    temp_paired_Pgt40_range = ScalingObj$temp_paired_Pgt40_range
    pquant = ScalingObj$pquant

    nboot <- length(ScalingBoot)-1  # first one is the raw unbootstrapped outcome

    for ( iboot in 1:nboot){

      ncb <- paste0("iboot",iboot)
      xtmp <- ScalingBoot[[iboot]]$xvar_mean
      ytmp <- ScalingBoot[[iboot]]$pr_quant

      if (iboot == 1){
        nx <- dim(ytmp)[1]
        nq <- dim(ytmp)[2]
        pr_quant_boot <- array(NA,dim= c(nx,nq,nboot))
        xvar_boot     <- array(NA,dim= c(nx,nboot))
      }
      pr_quant_boot[,,iboot] <- ytmp
      xvar_boot[,iboot] <- xtmp

    }

    pr_uncertainty <- apply(pr_quant_boot,FUN=quantile,probs=c(0.05,0.5,0.95),na.rm=TRUE,MARGIN = c(1,2))
    ltmp <- ! apply(! is.na(pr_quant_boot),FUN=all,MARGIN=c(1,2))

    for (ip in 1:3) {
      ytmp <-  pr_uncertainty[ip,,]
      ytmp[ltmp] <- NA
      pr_uncertainty[ip,,] <- ytmp
    }


    pcha <- c(16, 18,16, 15, 16,18,16 )
    cexa <- c(0.0, 1.2,1.5, 1.8, 1.5,1.2 ,0.0)

    yline <- 1.1*rep(ylimits[1], length(temp_paired_range))
    if (is == 2) {yline <- 1.2*rep(ylimits[1],length(temp_paired_range))}
    pyfac <- 1.25
    lines(temp_paired_wet_range,yline,col=lcols[1,is],lwd=2)
    points(temp_paired_wet_range,yline,pch = pcha,
           col = lcols[1,is],cex = 0.6*cexa)
    if (is == 1){ text(temp_paired_wet_range[7],yline[7],">0.1",pos=4,cex=0.7) }

    yline <- pyfac*yline
    lines(temp_paired_Pgt5_range,yline,col=lcols[2,is],lwd=2)
    points(temp_paired_Pgt5_range,yline,pch = pcha,
           col = lcols[2,is],cex =  0.6*cexa)
    if (is == 1){text(temp_paired_Pgt5_range[1],yline[1],">5",pos=2,cex=0.7) }

    yline <- pyfac*yline
    lines(temp_paired_Pgt20_range,yline,col=lcols[3,is],lwd=2)
    points(temp_paired_Pgt20_range,yline,pch = pcha,
           col = lcols[3,is],cex =  0.6*cexa)
    if (is == 1){ text(temp_paired_Pgt20_range[1],yline[1],">20",pos=2,cex=0.7) }


    # plot bootstrap results and fits

    for (ii in 1:3) {

      ltsel <- !is.na(pr_quant[,ii]) & xvar_mean >= temp_paired_Pgt5_range[1] & xvar_mean <= temp_paired_Pgt5_range[7]
      if (loessrange[1] != -999 & loessrange[2] - loessrange[1] > 5  )
      {ltsel <- !is.na(pr_quant[,ii]) & xvar_mean > loessrange[1] & xvar_mean <loessrange[2]}

      if (any(ltsel)) {
        if  (range(xvar_mean[ltsel])[2] - range(xvar_mean[ltsel])[1] > 5) {

          xx <-   c(xvar_mean,rev(xvar_mean))
          yy <- c(pr_uncertainty[1,,ii],rev(pr_uncertainty[3,,ii]))
          #  xx <-   c(xvar_mean[ltsel],rev(xvar_mean[ltsel]))
          #  yy <- c(pr_uncertainty[1,ltsel,ii],rev(pr_uncertainty[3,ltsel,ii]))

          xx <- xx[!is.na(yy)   ]
          yy <- yy[!is.na(yy)   ]
          polygon(xx,yy,col=t_col(ucols[is],percent=ucoltrans[is]),border=NA)
          rlspan = 0.5

            print("x interpolated fit")
            yyin <- pr_quant[ltsel,ii]
            xxin <- xvar_mean[ltsel]
               xxout <- seq(2,26,0.1)
            xfit <- gl.smooth(xin=xxin,yin=yyin,rlspan=rlspan, xout=xxout,drange=range(xxout))
            lines(xfit$x,xfit$y,col=pctcols[ii,is],lwd=3)

            if (lpoutsideonly) {

			pdiff <- (pr_uncertainty[3,,ii]-pr_uncertainty[1,,ii])/ pr_uncertainty[2,,ii]
			ltmp <- pdiff > 0.2
			ltmp[is.na(ltmp)] <- FALSE
		    ltmp = ! ltsel
		    ltmp <- ltmp | c(ltmp[-1],FALSE) | c(FALSE,ltmp[-length(ltmp)])
		    #  ltmp <- ltmp | c(ltmp[-1],FALSE) | c(FALSE,ltmp[-length(ltmp)])
		    ltmp <- ltmp | is.na(pdiff)


		      # expand a little

              points(xvar_mean[ltmp],pr_quant[ltmp,ii],col=pctcols[ii,is],pty="s",pch=14+ii)
			} else {
               points(xvar_mean,pr_quant[,ii],col=pctcols[ii,is],pty="s",pch=14+ii)
		    }

          } # if
        } # if

      } # over percentiles ii



    }

xfac = 0.08*(xlimits[2]-xlimits[1])
text(xlimits[1]+0.6*xfac,0.83*ylimits[2],plottitlelt[1],pos=4,cex=1.1)
text(xlimits[1]+0.6*xfac,0.62*ylimits[2],plottitlelt[2],pos=4,cex=1.1)
xp <- c(xlimits[1]+0.15*xfac,xlimits[1]+0.4*xfac, xlimits[1]+0.65*xfac)
points(xp,rep(0.83*ylimits[2],3),pch=c(15,16,17),col=pctcols[,1])
points(xp,rep(0.62*ylimits[2],3),pch=c(15,16,17),col=pctcols[,2])
text(xlimits[1]-1.5,0.83*ylimits[2],paste0(plabel,")"),pos=2,cex=1.1,xpd=NA)

dev.off()

}

## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  return(t.col)
}
## END



