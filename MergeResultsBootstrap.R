#
# note T is used for dew point, Tas for temperature

# master script to make all plots, replaces the nonbootstrapped version
#

# # custum settings

luncertainty <- TRUE
labspercentile <- FALSE
plotset <- 1  # 1 NL 3 med en 4 cent SFR
CPprname <- "pr_mean"
CPprname <- "pr_sample"
seasonid <- "MJJAS"


hdir <- Sys.getenv("HOME")
#source(paste0(hdir,"/mysoft/R_lib/mylib.R"))
#source(paste0(hdir,"/mysoft/R_lib/ScalingLib.R"))
wdir <- paste0(hdir,"/analysis/PRINCIPLES/analysis_github/")

source(paste0(wdir,"/ScalingLibSmall.R"))
iepsout = 1 # force eps file

setwd(wdir)

library("ncdf4",         lib.loc=Rlibdir)
library("ncdf4.helpers", lib.loc=Rlibdir)
library("abind",         lib.loc=Rlibdir)
library("rlist",         lib.loc=Rlibdir)


plotopt <- list(IntDiffmidRHplot = c(-60,130), 
                IntDifflowRHplot = c(-60,130),
                whfhigh = c(0,0.25),
                whfmid = c(0,0.2),
                whflow = c(0,0.15),
                whfall = c(0,0.25))

if ( labspercentile) {
  plotopt <- list(IntDiffmidRHplot = c(-90,60), 
                  IntDifflowRHplot = c(-90,60),
                  whfhigh = c(0,0.25),
                  whfmid = c(0,0.2),
                  whflow = c(0,0.1),
                  whfall = c(0,0.25))
}



for (plotset in c(1,3,4)) {    
    
    # four data sets 
    if ( plotset == 1) {
      areaid <- "NLstations.land.all"
      cobs   <- "OBS-NL-STATIONS-NLstations"
      cobs2   <- "OBS_NL"
    }
    
    if ( plotset == 2) {
      cobs <- "OBS-SFR-STATIONS-SFRstations"
      cobs2 <- "OBS_SFR"
      areaid <- "SFRstations.land.all"
    }
    
    if ( plotset == 3) {
      cobs <- "OBS-MEDSEA-STATIONS-MEDSEAstations"
      cobs2 <- "OBS_SFR-med"
      areaid <- "SFRstationsMedSea.land.all"
    }
    if ( plotset == 4) {
      cobs <- "OBS-NoMEDSEA-STATIONS-NoMEDSEAstations"
      cobs2 <- "OBS_SFR-cent"
      areaid <- "SFRstationsNoMedSea.land.all"
    }
    
    
    # directories in and output
    dirdata <-  "data_hourly_noboot/"
    dirdata_boot <-  "data_hourly_bootstrap/"
    dirplots <- "plots_merged_bootstrap/" 
    
    if ( ! dir.exists(dirplots)) { dir.create(dirplots)}
    
    evalp <- "1991-2010"
    pto <- 19
    ptc <- 5
    ptr <- 8
    
    EXPLIST = list(
      
          # for the moment three times the same data here 
          # OBS 
          # OBS   =    list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "black",  pt=pto, lw = 3,  pname = "OBS 1991-2020"),
          # OBS2   =   list(period= "2011-2020", modelid = cobs,       prname = "precip", col = "grey40",   pt=pto, lw = 3, pname = "OBS 2011-2020"),
          # OBS3   =   list(period= "1991-2010", modelid = cobs,       prname = "precip", col = "grey80",   pt=pto, lw = 3,  pname = "OBS 1991-2010"),
          OBS   =    list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "black",  pt=pto, lw = 3,  pname = "OBS 1991-2020"),
          OBS2   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey40",   pt=pto, lw = 3, pname = "OBS 1991-2020"),
          OBS3   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey80",   pt=pto, lw = 3,  pname = "OBS 1991-2020"),
          
          # CPMs
          ETH =      list(period= "CTL", modelid = "CPM-ERAINT-ETH-ALP3",         prname = CPprname, col = "red",    pt=ptc, lw = 1.5,  pname = "COSMO-CLM"),
          METO =     list(period= "CTL", modelid = "CPM-ERAINT-METO-ALP3",        prname = CPprname, col = "magenta",pt=ptc, lw =  1.5, pname = "UKMO-UM"),
          HCLIMALP3 =list(period= "CTL", modelid = "CPM-ERAINT-HCLIM-ALP3",       prname = CPprname, col = "blue",   pt=ptc, lw =  1.5, pname = "HCLIM-ALP"),
          HCLIMNWE3 =list(period= "CTL", modelid = "CPM-ERAINT-HCLIM-NWE3",       prname = CPprname, col = "cyan",   pt=ptc, lw =  1.5, pname = "HCLIM-NWE"),
          CNRM      =list(period= "CTL", modelid = "CPM-ERAINT-CNRM-NWE3",        prname = CPprname, col = "green",  pt=ptc, lw =  1.5, pname = "AROME41"),
          
          # CRCMs
          CLMCOM =   list(period= evalp, modelid = "RCM-ERAINT-CLMcom-ETH-COSMO", prname = "pr",     col = "red",   pt=ptr, lw = 2,    pname = "COSMO-CLM*"),
          HADRM =    list(period= evalp, modelid = "RCM-ERAINT-MOHC-HadREM3",     prname = "pr",     col = "magenta",   pt=ptr, lw = 2,pname = "HadRM3"),
          RACMO =    list(period= evalp, modelid = "RCM-ERAINT-KNMI-RACMO22E",    prname = "pr",     col = "blue",pt=ptr, lw = 2,      pname = "RACMO"),
          RCAO =     list(period= evalp, modelid = "RCM-ERAINT-SMHI-RCA4",        prname = "pr",     col = "cyan",   pt=ptr, lw = 2,   pname = "RCA4"),
          HIRHAM   = list(period= evalp, modelid = "RCM-ERAINT-DMI-HIRHAM5",      prname = "pr",     col = "green4" ,   pt=ptr, lw = 2,pname = "HIRHAM5"),
          REMO =     list(period= evalp, modelid = "RCM-ERAINT-GERICS-REMO2015",  prname = "pr",     col = "yellow2",   pt=ptr, lw = 2,pname = "REMO"),
          ALADIN =   list(period= evalp, modelid = "RCM-ERAINT-CNRM-ALADIN63",    prname = "pr",     col = "green",   pt=ptr, lw = 2,  pname = "ALADIN")
          
    )
    
    
    
    xpdf <- c(0.2,0.5,0.8)
    npdf <- length(xpdf)
    ipdf <- array(NA,npdf)
    
    nexps <- length(EXPLIST)
    iexps <-  1:nexps
    #iexps <-  13
    
    ccols <-c()
    pchs <- c()
    
    #===================================================================================================
    gl.FitScalingRange <- function(ScalingObj=ScalingObj,temprange = c(6,26),rlspan=0.6, labspercentile = FALSE, lintensitysel = TRUE){
      
      nbins = ScalingObj$nbins
      xvar_mean = ScalingObj$xvar_mean
      pr_quant = ScalingObj$pr_quant
      NFreqExc = ScalingObj$NFreqExc
      ExcTres  = ScalingObj$ExcTres
      pquant = ScalingObj$pquant
      temp_paired_Pgt5_range = ScalingObj$temp_paired_Pgt5_range
        
      nquant = length(pquant)
     # print(paste("nq =",nquant))
      
      zalpha_lm <- array(NA,nquant)
      zbeta_lm <- array(NA,nquant)
      Pmin_lm <- array(NA,nquant)
      Pmax_lm <- array(NA,nquant)
      Tmin_lm <- array(NA,nquant)
      Tmax_lm <- array(NA,nquant)
      FitErr_lm <- array(NA,nquant)
      SumErr_lm <- array(NA,nquant)
      AbsSumErr_lm <- array(NA,nquant)
      nSum_lm <- array(NA,nquant)
    
      tempLOESS <- seq(0,30,0.1)
      prLOESS  <-  array(NA,c(length(tempLOESS),nquant))
      
      # weights with fitting
      ltsel <-  xvar_mean > temp_paired_Pgt5_range[1] & xvar_mean <temp_paired_Pgt5_range[7]
      for (ip in 1:nquant){
        
          xvar_min <- temprange[1]
          xvar_max <- temprange[2]
          lsel <-  xvar_mean > xvar_min & xvar_mean < xvar_max & !is.na(pr_quant[,ip]) & pr_quant[,ip] > 0.1
          if (lintensitysel) lsel <- lsel & ltsel
            
          if ( sum(! is.na(pr_quant[lsel,ip])) >= 7 ) {
            
            LinMod <- lm(log(pr_quant[lsel,ip])~xvar_mean[lsel],na.action = na.omit)
            zalpha_lm[ip]   <- LinMod$coefficients[2]
            zbeta_lm[ip]    <- LinMod$coefficients[1]
            Pmin_lm[ip] <- exp(zbeta_lm[ip]+zalpha_lm[ip]*xvar_min)
            Pmax_lm[ip] <- exp(zbeta_lm[ip]+zalpha_lm[ip]*xvar_max)
            Tmin_lm[ip] <- xvar_min
            Tmax_lm[ip] <- xvar_max
            FitErr_lm[ip] <- sd(exp(LinMod$residuals)-1)
            SumErr_lm[ip] <- sum(exp(LinMod$residuals)-1)
            AbsSumErr_lm[ip] <- sum(abs(exp(LinMod$residuals)-1))
            nSum_lm[ip] <- length(LinMod$residuals)
          }
          
          # loesfilter
          if ( sum(! is.na(pr_quant[lsel,ip])) >= 7 ) {
     #       print(paste(ip,range(xvar_mean[lsel])))
            prLOESS[,ip] <- gl.smooth(xin=xvar_mean[lsel],yin=pr_quant[lsel,ip],rlspan=rlspan, xout=tempLOESS,drange=range(xvar_mean[lsel]))$y
          }
      }
      
      # some addional variables
      
      lsel <- ScalingObj$Nall > 1000 & ltsel
      frac_wet <- ScalingObj$Nwet/ScalingObj$Nall
      # print(frac_wet)
      rspanfrac = 0.3
      if (labspercentile) {
        frac_wet_smooth <- gl.smooth(xin=xvar_mean[lsel],yin=frac_wet[lsel],rlspan= rspanfrac, xout=tempLOESS,drange=range(xvar_mean[lsel]))$y
      } else {
        frac_wet_smooth <- gl.smooth(xin=xvar_mean[lsel],yin=frac_wet[lsel],rlspan= rspanfrac, xout=tempLOESS,drange=range(xvar_mean[lsel]))$y
      }
      
      
      nall <- ScalingObj$Nall
      
      
      FitErr <- list("StErr" = FitErr_lm, "SumErr" = SumErr_lm, "AbsErr" = AbsSumErr_lm, "npfit" = nSum_lm)
      FitObj <- list("FitErr" = FitErr )
      FitObj$zalpha_lm = exp(zalpha_lm)-1
      FitObj$zbeta_lm  = zbeta_lm
      FitObj$Pmin_lm  <- Pmin_lm
      FitObj$Pmax_lm  <- Pmax_lm
      FitObj$Tmin_lm  <- Tmin_lm
      FitObj$Tmax_lm  <- Tmax_lm
      FitObj$xLOESS   <- tempLOESS
      FitObj$prLOESS  <- prLOESS
      FitObj$fwet     <- frac_wet_smooth
      FitObj$fwet_noint   <- frac_wet
      FitObj$temp_noint   <- xvar_mean
      FitObj$nall_noint   <- nall
      FitObj$origscalingobj <- ScalingObj
    
      return(FitObj)
      
    }
    
    
    
    # read in data in a number of lists
    
    
    FitHighRH <- list()
    FitMidRH<- list()
    FitLowRH <- list()
    FitAll <- list()
    FitTasAll <- list()
    ExtremesDataList <- list()
    CDFDataList <- list()
    
    
    for (iexp in iexps){
      
        ccols <- c(ccols,EXPLIST[[iexp]]$col)
        pchs <-  c(pchs, EXPLIST[[iexp]]$pt)
        varnameP <- EXPLIST[[iexp]]$prname
        
        if (iexp %in% 1:3){ # obs data set
          
          ddid = paste(varnameP,cobs,seasonid,"1991-2020",sep="-")
          
        } else {
          
          ddid = paste(varnameP,EXPLIST[[iexp]]$modelid,areaid,seasonid,EXPLIST[[iexp]]$period,sep="-")
          
        }
        
    
        dataid <- "Cluster0_RHsens"
        Rdatafile <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
        
        dataid <- "Cluster"
        Rdatafile2 <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
        
         
        
        print(Rdatafile)
        
        load(Rdatafile)
        load(Rdatafile2)
        
        rspan = 0.6
        
        FitData <- gl.FitScalingRange(ScalingData_high,rlspan=rspan)
        FitHighRH  <- list.append(FitHighRH,FitData)
        
        FitData <- gl.FitScalingRange(ScalingData_mid,rlspan=rspan)
        FitMidRH  <- list.append(FitMidRH,FitData)
        FitData <- gl.FitScalingRange(ScalingData_low,rlspan=rspan)
        FitLowRH  <- list.append(FitLowRH,FitData)
        FitData <- gl.FitScalingRange(ScalingList[["all"]],rlspan=rspan)
        FitAll <- list.append(FitAll,FitData)
        
        FitData <- gl.FitScalingRange(ScalingListTas[["all"]],temprange = c(10,30),rlspan=rspan)
        FitTasAll <- list.append(FitTasAll,FitData)
        
        ExtremesDataList <-list.append(ExtremesDataList,ExtremesList[["all"]])
        CDFDataList <-list.append(CDFDataList,CDFlist[["all"]])
        
    }
    
    # xx <- ScalingList$all$xvar_mean
    # yy <- ScalingList$all$pr_quant[,2]
    # 
    # plot(xx,yy,log="y",xlim=c(4,26))
    # lines(FitAll[[1]]$xLOESS,FitAll[[1]]$prLOESS[,2])
    
    ResultsBootstrap <- new.env()
    
    if (luncertainty) {
      
          varnameP <- "precip"  # OBS-NL-STATIONS-NLstations-SHY-1991-2020.Rdata*
          dataid = paste(varnameP,cobs,seasonid,"1991-2020",sep="-")
          fileobj <- paste0(dirdata_boot,"/Cluster0_AllStatBootstrap_",dataid,".Rdata")
          load(file=fileobj, envir = ResultsBootstrap)
          
      
          # compute some statistics from bootstrapped observations
          
          ScalingListObs    <- ResultsBootstrap$ScalingList
          ScalingListObsFit <- list()
          
          nboot <- length(ScalingListObs) - 1
          for ( iboot in 1:(nboot+1)){
            ctmp <- names(ScalingListObs)[iboot]
            ScalingListObsFit[[ctmp]] <- gl.FitScalingRange(ScalingListObs[[iboot]],rlspan=rspan)
            ScalingListObsFit[[ctmp]]$RHhigh <- gl.FitScalingRange(ScalingListObs[[iboot]]$RHhigh,rlspan=rspan)
            ScalingListObsFit[[ctmp]]$RHlow <- gl.FitScalingRange(ScalingListObs[[iboot]]$RHlow,rlspan=rspan)
            ScalingListObsFit[[ctmp]]$RHmid <- gl.FitScalingRange(ScalingListObs[[iboot]]$RHmid,rlspan=rspan)
            # CGL
            ScalingListObsFit[[ctmp]]$RHall  <- gl.FitScalingRange(ScalingListObs[[iboot]],rlspan=rspan)
            
          }
    
    }
    
    
    #===================================================================================================
    gl.UncertaintyBoot <- function(ListIn,StatName,nboot=100,LevelName = NA, LevelName2 = NA, rtype = "rel",nbootmargin=10) {
      
      # ListIn <- ScalingListObsFit
      # StatName <- "prLOESS"
      # nboot=100
      # LevelName = NA
      # LevelName2 = NA
      
      ctmp = "raw"
      
      if (! is.na(LevelName)) {
        Xstat <- ListIn[[ctmp]][[LevelName]][[StatName]]
        if (! is.na(LevelName2)) {
          if (rtype == "rel") { Xstat <- 100*(Xstat /  ListIn[[ctmp]][[LevelName2]][[StatName]] - 1)} 
          if (rtype == "abs") { Xstat <- Xstat -  ListIn[[ctmp]][[LevelName2]][[StatName]] } 
        }
      } else {
        Xstat <- ListIn[[ctmp]][[StatName]]
      }
      
      ndim <- dim(Xstat)
      if (is.null(ndim)) {ndim=length(Xstat)}
      Xstat_boot <- array(NA,c(ndim,nboot))
      XstatNoBoot <-Xstat 
      
      for ( iboot in 1:nboot){
        ctmp = paste0("iboot",iboot)
        if (! is.na(LevelName)) {
          Xstat <- ListIn[[ctmp]][[LevelName]][[StatName]]
          if (! is.na(LevelName2)) {
            if (rtype == "rel") { Xstat <- 100*(Xstat /  ListIn[[ctmp]][[LevelName2]][[StatName]] - 1)} 
            if (rtype == "abs") { Xstat <- Xstat -  ListIn[[ctmp]][[LevelName2]][[StatName]] } 
          }
        } else {
          Xstat <- ListIn[[ctmp]][[StatName]]
        }
        if (length(ndim) == 1){Xstat_boot[,iboot]   <- Xstat}
        if (length(ndim) == 2){Xstat_boot[,,iboot]  <- Xstat}
        if (length(ndim) == 3){Xstat_boot[,,,iboot] <- Xstat}
      }
      
      xtmp <- apply(Xstat_boot,FUN=quantile,probs=c(0.05,0.5,0.95),na.rm=TRUE,MARGIN = 1:length(ndim))   
      ltmp <- ! apply(! is.na(Xstat_boot),FUN=all,MARGIN=1:length(ndim))
      
      if ( ! is.na(nbootmargin)) {
          lmask <- ! is.na(Xstat_boot) 
          lsum <- apply(lmask,FUN=sum,MARGIN=1:length(ndim))
          ltmp <- lsum < nboot-nbootmargin
      }
      
      if (length(ndim) == 1){
      for (ip in 1:3) {
        ytmp <-  xtmp[ip,]
        ytmp[ltmp] <- NA
        xtmp[ip,] <- ytmp
      }
      XstatP05 <- xtmp[1,]  
      XstatP50 <- xtmp[2,]  
      XstatP95 <- xtmp[3,]  
      }
      
      if (length(ndim) == 2){
        for (ip in 1:3) {
          ytmp <-  xtmp[ip,,]
          ytmp[ltmp] <- NA
          xtmp[ip,,] <- ytmp
        }
        XstatP05 <- xtmp[1,,]  
        XstatP50 <- xtmp[2,,]  
        XstatP95 <- xtmp[3,,]  
      }
    
      if (length(ndim) == 3){
        for (ip in 1:3) {
          ytmp <-  xtmp[ip,,,]
          ytmp[ltmp] <- NA
          xtmp[ip,,,] <- ytmp
        }
        XstatP05 <- xtmp[1,,,]  
        XstatP50 <- xtmp[2,,,]  
        XstatP95 <- xtmp[3,,,]  
      }
      
      BootRange <- list()
      BootRange[[StatName]]$NoBoot <- XstatNoBoot
      BootRange[[StatName]]$P05 <- XstatP05
      BootRange[[StatName]]$P50 <- XstatP50
      BootRange[[StatName]]$P95 <- XstatP95
      return(BootRange)
      
      }
    
    
    if (luncertainty){
        yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = "prLOESS")
        xx <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = "xLOESS")
        FitAll$uncertainty$xLOESS <- xx$xLOESS
        FitAll$uncertainty$prLOESS <- yy$prLOESS
    } 
    
    #===================================================================================================
    gl.PlotScalingComparison <- function(ScalingList,ScalingList2 = NA,explist,psel=0,
                                         plotfile = "dummy.eps",ptype=2, iquant=3, ytitle =  "no", 
                                         ilegend = 0, ptext = "", xlimits = c(4,23), ylimforced = NA,
                                         xlabel="dew point temperature [degC]", luncertainty = FALSE) {
      
      
      ip <- iquant
      if (psel[1]==0){psel <- 1:length(explist)}
    
      fileout <- plotfile
      source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
      options(scipen=2)
      par(cex=1.4)
      
      
      hlines = c(0.5,1,2,5,10,20,50)
      vlines <- c(5,10,15,20,25)
      
      lcc <- FALSE
      
      if ( ptype == 1){
        
        cyplot <- "ScalingList[[isel]]$prLOESS[,ip]"
        cxplot <- "ScalingList[[isel]]$xLOESS"
        cy2plot1 <- "ScalingList[['uncertainty']]$prLOESS$P05[,ip]"
        cy2plot2 <- "ScalingList[['uncertainty']]$prLOESS$P95[,ip]"
        
        ylabel <- "Intensity [mm/hour]"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits <- c(1,60)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
        
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab=xlabel,log = "y")
    
        if (luncertainty){
          
          isel <- 1
          xx <- eval(parse(text=cxplot))
          yy1 <- eval(parse(text=cy2plot1))
          yy2 <- eval(parse(text=cy2plot2))
          x <-   c(xx,rev(xx))
          y <- c(yy1,rev(yy2))
          x <- x[!is.na(y)   ]
          y <- y[!is.na(y)   ]
          polygon(x,y,col="grey88",border=NA)
         # print(yy2-yy1)
        }
        
        lcc <- TRUE
        
      } 
      if ( ptype == 2){
        
        cyplot <- "100*(ScalingList2[[isel]]$prLOESS[,ip]/ScalingList[[isel]]$prLOESS[,ip]-1)"
        cxplot <- "ScalingList[[isel]]$xLOESS"
        ylabel = "Intensity difference 1 vs 2 [%]"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits = c(-40,60)
        hlines <- c(-40,-20,0,20,40,60,80,100)
        hlines <- c(-50,0,50,100)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
        
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab=xlabel)
        for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
        for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
        
        if (luncertainty){
          
          yy_1_low <- ScalingList[['uncertainty']]$prLOESS$P05[,ip]
          yy_1_up  <- ScalingList[['uncertainty']]$prLOESS$P95[,ip]
          yy_1_mid  <- ScalingList[['uncertainty']]$prLOESS$P50[,ip]
          yy_2_low <- ScalingList2[['uncertainty']]$prLOESS$P05[,ip]
          yy_2_mid  <- ScalingList2[['uncertainty']]$prLOESS$P50[,ip]
          yy_2_up  <- ScalingList2[['uncertainty']]$prLOESS$P95[,ip]
          
          isel <- 1
          xx <- eval(parse(text=cxplot))
          yy1 <- 100.*(yy_2_up / yy_1_low  - 1)
          yy2 <- 100.*(yy_2_low / yy_1_up  - 1)
          yym <- 100.*(yy_2_mid/ yy_1_mid - 1)
          yy1 <- yym + (yy1-yym)/1.4
          yy2 <- yym + (yy2-yym)/1.4
          
          x <-   c(xx,rev(xx))
          y <- c(yy1,rev(yy2))
          x <- x[!is.na(y)   ]
          y <- y[!is.na(y)   ]
          polygon(x,y,col="grey88",border=NA)
        }
        
        
      } 
      if ( ptype == 3){
        cyplot <- "ScalingList[[isel]]$prLOESS[,ip]"
        cy2plot <- "ScalingList2[[isel]]$prLOESS[,ip]"
        cxplot <- "ScalingList[[isel]]$xLOESS"
        ylabel <- "Intensity [mm/hour]"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits <- c(1,60)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
        
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab=xlabel,log = "y")
      } 
      
      if ( ptype == 4){
        
        cyplot <- "ScalingList[[isel]]$prLOESS[,ip]"
        cxplot <- "ScalingList[[isel]]$xLOESS"
        ylabel <- "Intensity [mm/hour]"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits <- c(1,60)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
        lcc <- TRUE
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab=xlabel,log = "y")
        
        
      } 
      
      if ( ptype == 5){
        
        cyplot <- "ScalingList2[[isel]]$fwet/ScalingList[[isel]]$fwet"
        cxplot <- "ScalingList[[isel]]$xLOESS"
        ylabel <- "fwet high humidity versus low"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits <- c(0,2)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
     #   xlimits = c(4,24)
        
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab=xlabel)
        
      } 
      
     
      if ( ptype == 6){
        
        cyplot <- "ScalingList[[isel]]$fwet"
        cxplot <- "ScalingList[[isel]]$xLOESS"
        cy2plot <- "ScalingList[['uncertainty']]$fwet"
        
        ylabel <- "fwet"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits <- c(0,0.25)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
        #xlimits = c(4,24)
        
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab=xlabel)
        
        if (luncertainty){
          isel <- 1
          xx <- eval(parse(text=cxplot))
          yy <- eval(parse(text=cy2plot))
          x <-   c(xx,rev(xx))
          y <- c(yy$P05,rev(yy$P95))
          x <- x[!is.na(y)   ]
          y <- y[!is.na(y)   ]
          polygon(x,y,col="grey88",border=NA)
        }
        
        
        grid()
        
      } 
      
      
      # actual plotting
      
      pcol <- c()
      plegend <- c()
      
      if (lcc){
        CCrate <- 0.065
        xcommon <- seq(-5,30,1)
        yscaling7 <- 2*(1+CCrate)**xcommon # plot some reference lines CC and 2CC
        yscaling14 <- (1+2*CCrate)**xcommon
        lines(xcommon,yscaling14,lty=3,lwd=1,col="red")
        lines(xcommon,2*yscaling14,lty=3,lwd=1,col="red")
        lines(xcommon,yscaling7,lty=3,lwd=1,col="black")
        lines(xcommon,2*yscaling7,lty=3,lwd=1,col="black")
      }
      
      for (isel in psel) {
        xx <- eval(parse(text=cxplot))
        yy <- eval(parse(text=cyplot))
        
        pcol <- c(pcol,EXPLIST[[isel]]$col)
        plegend <- c(plegend,EXPLIST[[isel]]$pname)
        
        lines(xx,yy,col=EXPLIST[[isel]]$col,lw=EXPLIST[[isel]]$lw)
        if (ptype == 3){
          yy2 <- eval(parse(text=cy2plot))
          lines(xx,yy2,col=EXPLIST[[isel]]$col,lw=EXPLIST[[isel]]$lw,lty=2)
        }
      }
      
      if ( ptext != "" ){
        if (ilegend == 0) { text(xlimits[2],ylimits[1],ptext,pos=2 )}
        if (ilegend == 1) { text(xlimits[1],ylimits[2],ptext,pos=4 )}
        if (ilegend == 2) { text(xlimits[2],ylimits[1],ptext,pos=2)}
      }
      
      if  ( ilegend != 0 ) {
        if (ilegend == 1) {legend(xlimits[2],ylimits[1],plegend,col=pcol,xjust=1,yjust=0,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
        if (ilegend == 2) {legend(xlimits[1],ylimits[2],plegend,col=pcol,xjust=0,yjust=1,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
      }
      
     
      dev.off()
      
    }
    
    
    #===================================================================================================
    gl.PlotCDFextremesComparison <- function(ExtremesList,explist,psel=0,
                                         plotfile = "dummy.eps",ptype=1, ytitle =  "no", ilegend = 0, ptext = "", xlimits = c(1,5e-8), 
                                         ylimforced = NA, luncertainty = FALSE) {
      
      if (psel[1]==0){psel <- 1:length(explist)}
      fileout <- plotfile
      source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
      options(scipen=2)
      par(cex=1.4)
      
      
      hlines = c(10,20,30,50,70)
      vlines <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)
      
      if ( ptype == 1){
        cyplot <- "ExtremesList[[isel]]$pvar_extremes"
        cy2plot <- "ExtremesList[['uncertainty']]$pvar_extremes"
        cxplot <- "1-ExtremesList[[isel]]$pctls_extremes"
        ylabel <- "Intensity [mm/hour]"
        if (ytitle != "no")  { ylabel = ytitle }
        ylimits <- c(1,80)
        if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
        
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="probability of exceedance",log = "x")
        for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
        for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
        
        
      } 
      pcol <- c()
      plegend <- c()
      
       isel <- psel[1]
       xx <- eval(parse(text=cxplot))
       
       if (luncertainty){
           yy <- eval(parse(text=cy2plot))
           x <-   c(xx,rev(xx))
           y <- c(yy$P05,rev(yy$P95))
           
           y [x>0.04] <- NA
           x <- x[!is.na(y)   ]
           y <- y[!is.na(y)   ]
           polygon(x,y,col="grey88",border=NA)
       }
       
      for (isel in psel) {
        
        xx <- eval(parse(text=cxplot))
        yy <- eval(parse(text=cyplot))
        
        pcol <- c(pcol,EXPLIST[[isel]]$col)
        plegend <- c(plegend,EXPLIST[[isel]]$pname)
    
        lines(xx,yy,col=EXPLIST[[isel]]$col,lw=EXPLIST[[isel]]$lw)
        #  if (ptype == 3){
        #   yy2 <- eval(parse(text=cy2plot))
        #   lines(xx,yy2,col=EXPLIST[[isel]]$col,lw=EXPLIST[[isel]]$lw)
        # }
      }
      
      if ( ptext != "" ){
        if (ilegend == 0) { text(xlimits[2],ylimits[1],ptext,pos=2 )}
        if (ilegend == 1) { text(xlimits[1],ylimits[2],ptext,pos=4 )}
        if (ilegend == 2) { text(xlimits[2],ylimits[1],ptext,pos=2)}
      }
      
      if  ( ilegend != 0 ) {
        if (ilegend == 1) {legend(xlimits[2],ylimits[1],plegend,col=pcol,xjust=1,yjust=0,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
        if (ilegend == 2) {legend(xlimits[1],ylimits[2],plegend,col=pcol,xjust=0,yjust=1,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
      }
      
      
      dev.off()
      
    }
    
    
     
    # defines sets to be plotted
     
    def_olist <- list(
      "CPMs" = list("psel" = 1:8,
                   "cdum" =  paste0("CPM_",CPprname)),
      "RCMs" = list("psel" = c(1,9:15),
                   "cdum" =  "RCM_pr")
      )
     
    if (luncertainty){   # only plots one obs period
        def_olist <- list(
          "CPMs" = list("psel" = c(1,4:8),
                        "cdum" =  paste0("CPM_",CPprname)),
          "RCMs" = list("psel" = c(1,9:15),
                        "cdum" =  "RCM_pr")
        )
    }
    
    
    for (io in 1:2) {
      
     cdum <- def_olist[[io]]$cdum  
     psel <- def_olist[[io]]$psel
     cmodel <- names(def_olist)[io]
     cdum <- paste0(cdum,ldash,seasonid)
     
     # plot not depedendent on quantile
     pptext = paste(cmodel,strsplit(cobs2,"_")[[1]][2],seasonid)
     
     # plot probability of exceedance plot
     
     # add uncertainty to list
     if (luncertainty){
         xvar <- "pvar_extremes"
         ExtremesListObs    <- ResultsBootstrap$ExtremesList
         yy <- gl.UncertaintyBoot(ListIn=ExtremesListObs,StatName = xvar)
         ExtremesDataList$uncertainty[[xvar]] <- yy[[xvar]]
     } 
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_CDF.eps")
     gl.PlotCDFextremesComparison(ExtremesList=ExtremesDataList,psel=psel,
                                  plotfile = plotfile,ilegend=2,ptext = pptext, luncertainty = luncertainty, xlimits = c(1,5e-7))
     
     
     # add uncertainty to list
     if (luncertainty){
         xvar <- "fwet"
         yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHhigh")
         FitHighRH$uncertainty[[xvar]] <- yy[[xvar]]
         yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHlow")
         FitLowRH$uncertainty[[xvar]] <- yy[[xvar]]
         yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHmid")
         FitMidRH$uncertainty[[xvar]] <- yy[[xvar]]
         # CGL
         yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHall")
         FitAll$uncertainty[[xvar]] <- yy[[xvar]]
         
     }
     
     
     ptype = 5 # difference in whf between low and high RH
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_ptype",ptype,".eps")
     gl.PlotScalingComparison(FitHighRH,FitLowRH,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype)
     
     ptype = 6
    
     pptext = paste(cmodel,strsplit(cobs2,"_")[[1]][2])
     
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_ptype",ptype,"_whf_highRH.eps")
     gl.PlotScalingComparison(FitHighRH,NA,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,
                              ytitle = "whf @ high RH", ylimforced = plotopt$whfhigh,ptext = pptext, ilegend = 0, luncertainty = luncertainty)
     
     
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_ptype",ptype,"_whf_lowRH.eps")
     gl.PlotScalingComparison(FitLowRH,NA,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,
                              ytitle = "whf @ low RH",ylimforced = plotopt$whflow,ptext = pptext, ilegend = 2, luncertainty = luncertainty)
     
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_ptype",ptype,"_whf_midRH.eps")
     gl.PlotScalingComparison(FitMidRH,NA,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,
                              ytitle = "whf @ mid RH",ylimforced = plotopt$whfmid,ptext = pptext, ilegend = 2, luncertainty = luncertainty)
     
    
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_ptype",ptype,"_whf_all.eps")
     gl.PlotScalingComparison(FitAll,NA,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,
                              ytitle = "whf",ylimforced = plotopt$whfall,ptext = pptext, ilegend = 2, luncertainty = luncertainty)
     
     
     plotfile <- paste0(dirplots,cdum,dash,areaid,"_ptype",ptype,"_whf_tas.eps")
     gl.PlotScalingComparison(FitTasAll,NA,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,
                              ytitle = "wet hour fraction [0-1] ",ylimforced = plotopt$whfall,ptext = pptext, ilegend = 2,
                              xlimits = c(10,30),xlabel = "temperature [degC]")
     
     for (iq in 2:3) {
       
       cpquant <- paste0("P",100*ScalingData_high$pquant[iq]) # take from high RH, assuming all percentiles are the same (should be)
       
         
       ptype = 1  # does plotting of single percentile 
       
       if (luncertainty){
           xvar <- "prLOESS"
           yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHhigh")
           FitHighRH$uncertainty[[xvar]] <- yy[[xvar]]
           yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHlow")
           FitLowRH$uncertainty[[xvar]] <- yy[[xvar]]
           yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHmid")
           FitMidRH$uncertainty[[xvar]] <- yy[[xvar]]
       } 
       
       # xvar <- "prLOESS"  # some testing on differences, between uncertainty from difference, and uncertainty sep. / sqrt(2). See evernote, almost no difference
       # yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar,LevelName = "RHlow", LevelName2 = "RHhigh")
       # xx <- FitAll[[1]]$xLOESS
       # plot(xx,yy$prLOESS$NoBoot[,2],xlim=c(5,24),ylim=c(-40,120),type="l", col="red")
       # lines(xx,yy$prLOESS$P05[,2])
       # lines(xx,yy$prLOESS$P95[,2])
       
       pptext = paste0(cmodel," ",cpquant," high RH"," ",strsplit(cobs2,"_")[[1]][2])
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_RHhigh.eps")
       gl.PlotScalingComparison(FitHighRH,NA,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ilegend = 0,ptext = pptext, luncertainty = luncertainty ) #ilegend = 2
       
       pptext = paste0(cmodel," ",cpquant," mid RH"," ",strsplit(cobs2,"_")[[1]][2])
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_RHmid.eps")
       gl.PlotScalingComparison(FitMidRH,NA,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ilegend = 0,ptext = pptext , luncertainty = luncertainty )
       
       pptext = paste0(cmodel," ",cpquant," low RH"," ",strsplit(cobs2,"_")[[1]][2])
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_RHlow.eps")
       gl.PlotScalingComparison(FitLowRH,NA,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ilegend = 0,ptext = pptext, luncertainty = luncertainty  )
     
       if (luncertainty){
           xvar <- "prLOESS"
           yy <- gl.UncertaintyBoot(ListIn=ScalingListObsFit,StatName = xvar)
           FitAll$uncertainty[[xvar]] <- yy[[xvar]]
       }
       
       pptext = paste0(cmodel," ",cpquant," all RH"," ",strsplit(cobs2,"_")[[1]][2])
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_RHall.eps")
       gl.PlotScalingComparison(FitAll,NA,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ilegend = 1,ptext = pptext, luncertainty = luncertainty)
       
       # temperature depdendency
       
       pptext = paste0(cmodel," ",cpquant," all RH"," ",strsplit(cobs2,"_")[[1]][2])
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_tas.eps")
       gl.PlotScalingComparison(FitTasAll,FitTasAll,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,
                                ptype=ptype,ilegend = 0,ptext = pptext,xlimits = c(10,30),xlabel = "temperature [degC]")
       
       
       ptype = 2 # plotting of difference between low/middle and high RH
       
       # uncertainty here is based on independence between errors in low and high humidity (rescaled by 1.4)
       # can be done better by looking at the samples themselves, but needs some coding
       # TODO
    
       pptext = paste0(cmodel," ",cpquant," ",strsplit(cobs2,"_")[[1]][2])
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_RHmid_vs_high.eps")
       gl.PlotScalingComparison(FitHighRH,FitMidRH,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ilegend = 0,
                                ytitle = "intensity difference mid RH vs high RH [%]",ptext = pptext,ylimforced = plotopt$IntDiffmidRHplot, luncertainty = luncertainty)
    
       plotfile <- paste0(dirplots,cdum,ldash,cpquant,dash,areaid,"_ptype",ptype,"_RHlow_vs_high.eps")
       gl.PlotScalingComparison(FitHighRH,FitLowRH,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ilegend = 0,
                                ytitle = "intensity difference low RH vs high RH [%]",ptext = pptext,ylimforced = plotopt$IntDifflowRHplot, luncertainty = luncertainty)
       
       # ptype = 3  # plots high and mid RH in one plot; uggly
       # plotfile <- paste0(dirplots,"CPM_",CPprname,"_P",100*ScalingData_high$pquant[iq],dash,areaid,"_ptype",ptype,"_comp.eps")
       # gl.PlotScalingComparison(FitHighRH,FitMidRH,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype,ytitle = "intensity difference low RH and high RH [%]"
       # ptype = 4 # plots without grid
       # plotfile <- paste0(dirplots,"CPM_",CPprname,"_P",100*ScalingData_high$pquant[iq],dash,areaid,"_ptype",ptype,".eps")
       # gl.PlotScalingComparison(FitHighRH,FitMidRH,EXPLIST,psel=psel,plotfile = plotfile,iquant=iq,ptype=ptype)
    
    
     
     }  
      
    }  
     
    
     
    # try to use the CDFs of dewpoint depression, etc
     
     
     
     #===================================================================================================
     gl.PlotCDFdist <- function(CDFList,explist,psel=0,
                                       plotfile = "dummy.eps",ptype=1, ytitle =  "no", ilegend = 0, ptext = "", xlimits = c(1,9), 
                                ylimforced = NA, pctl=0.5, luncertainty = FALSE) {
       
       if (psel[1]==0){psel <- 1:length(explist)}
       fileout <- plotfile
       source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
       options(scipen=2)
       par(cex=1.4)
       
       
       hlines = 1:30
       vlines <- c(1,2,3,5,7)
    
    #  some data processing
       
       
       if ( ptype == 1){
         
         ipctl <- which.min((CDFList[[1]]$xseq-pctl)**2)
         cyplot <- "CDFList[[isel]]$tddepresDist[,ipctl]"
         cxplot <- "1:length(CDFList[[isel]]$prthres)"
         ylabel <- "Dew point depression [degC]"
         if (ytitle != "no")  { ylabel = ytitle }
         ylimits <- c(0,10)
         if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
         plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="rainfall intensity",xaxt="n")
         xx <- vlines
         axis(side=1,at=xx,labels=names(CDFList[[1]]$prthres[xx]))
         
         for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
         for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
         
         if (luncertainty){
           isel <- 1
           xx <- eval(parse(text=cxplot))
           yy2 <- CDFList[["uncertainty"]]$tddepresDist$P95[,ipctl] 
           yy1 <- CDFList[["uncertainty"]]$tddepresDist$P05[,ipctl]
           ybias <- CDFList[["uncertainty"]]$tddepresDist$NoBoot[,ipctl]- CDFList[[1]]$tddepresDist[,ipctl]
           # cgl 
           ybias = 0.
           x <-   c(xx,rev(xx))
           y <- c(yy2,rev(yy1)) - ybias
           x <- x[!is.na(y)   ]
           y <- y[!is.na(y)   ]
           polygon(x,y,col="grey88",border=NA)
           # lines(xx,ybias,col="brown",lwd=2)
           
         }
         
         
       } 
     
       if ( ptype == 2){
         
         ipctl <- which.min((CDFList[[1]]$xseq-pctl)**2)
         cyplot <- "CDFList[[isel]]$tdDist[,ipctl]"
         cxplot <- "1:length(CDFList[[isel]]$prthres)"
         ylabel <- "Dew point [degC]"
         if (ytitle != "no")  { ylabel = ytitle }
         ylimits <- c(5,20)
         if ( ! is.na(ylimforced[1]) ) {ylimits <- ylimforced }
         plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="rainfall intensity",xaxt="n")
         xx <- vlines
         axis(side=1,at=xx,labels=names(CDFList[[1]]$prthres[xx]))
         
         for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
         for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
         
         if (luncertainty){
           isel <- 1
           xx <- eval(parse(text=cxplot))
           yy2 <- CDFList[["uncertainty"]]$tdDist$P95[,ipctl]
           yy1 <- CDFList[["uncertainty"]]$tdDist$P05[,ipctl]
           ybias <- CDFList[["uncertainty"]]$tdDist$NoBoot[,ipctl]- CDFList[[1]]$tdDist[,ipctl]
           # cgl
           ybias <- 0.
           x <-   c(xx,rev(xx))
           y <- c(yy2,rev(yy1))-ybias
           x <- x[!is.na(y)   ]
           y <- y[!is.na(y)   ]
           polygon(x,y,col="grey88",border=NA)
          # lines(xx,ybias,col="brown",lwd=2) # checked and marginal differences
         }
         
         
       } 
       
      
       pcol <- c()
       plegend <- c()
       
       for (isel in psel) {
         
         xx <- eval(parse(text=cxplot))
         yy <- eval(parse(text=cyplot))
         
         pcol <- c(pcol,EXPLIST[[isel]]$col)
         plegend <- c(plegend,EXPLIST[[isel]]$pname)
         
         lines(xx,yy,col=EXPLIST[[isel]]$col,lw=EXPLIST[[isel]]$lw)
       }
       
       if ( ptext != "" ){
         if (ilegend == 0) { text(xlimits[2],ylimits[1],ptext,pos=2 )}
         if (ilegend == 1) { text(xlimits[1],ylimits[2],ptext,pos=4 )}
         if (ilegend == 2) { text(xlimits[2],ylimits[1],ptext,pos=2)}
       }
       
       if  ( ilegend != 0 ) {
         if (ilegend == 1) {legend(xlimits[2],ylimits[1],plegend,col=pcol,xjust=1,yjust=0,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
         if (ilegend == 2) {legend(xlimits[1],ylimits[2],plegend,col=pcol,xjust=0,yjust=1,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
       }
       
       
       dev.off()
       
     }
     
    
    
    
    #===================================================================================================
    gl.PlotCDFdistFullsingle <- function(CDFList,explist,psel=0,
                                plotfile = "dummy.eps",ptype=1, ytitle =  "no", 
                                ilegend = 0, ptext = "", ylimforced = c(-5,10), cstat = "all",
                                luncertainty = FALSE) {
      
      if (psel[1]==0){psel <- 1:length(explist)}
      fileout <- plotfile
      source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
      options(scipen=2)
      par(cex=1.4)
      
      
      hlines = seq(ylimforced[1],ylimforced[2],(ylimforced[2]-ylimforced[1])/5    )
      vlines <- seq(0,1,0.25)
      
      #  some data processing
      
      iplot <- which(names(CDFList[[1]]$prthres) == cstat)
      
      
      if ( ptype == 1){
        
        cyplot <- "CDFList[[isel]]$tddepresDist[iplot,]"
        cxplot <- "CDFList[[isel]]$xseq"
        ylabel <- "DPD"
        if (ytitle != "no")  { ylabel = cstat }
        ylimits <- ylimforced 
        xlimits <- c(0,1)
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="cum. prob [0:1]")
       
        for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
        for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
      } 
      
      if ( ptype == 2){
        
        cyplot <- "CDFList[[isel]]$tddepresDist[iplot,]-CDFList[[1]]$tddepresDist[iplot,]"
        cxplot <- "CDFList[[isel]]$xseq"
        ylabel <- "Anom. DPD wrt obs 1991-2020 [degC]"
        if (ytitle != "no")  { ylabel = cstat }
        ylimits <- ylimforced 
        xlimits <- c(0,1)
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="cum. prob [0:1]")
        
        grid()
        abline(h=0,col="grey80",lty=2)
        
        if (luncertainty){
          isel <- 1
          xx <- eval(parse(text=cxplot))
          yy2 <- CDFList[["uncertainty"]]$tddepresDist$P95[iplot,] - CDFList[[1]]$tddepresDist[iplot,]
          yy1 <- CDFList[["uncertainty"]]$tddepresDist$P05[iplot,] - CDFList[[1]]$tddepresDist[iplot,]
          yy2 <- CDFList[["uncertainty"]]$tddepresDist$P95[iplot,] - CDFList[["uncertainty"]]$tddepresDist$NoBoot[iplot,]
          yy1 <- CDFList[["uncertainty"]]$tddepresDist$P05[iplot,] - CDFList[["uncertainty"]]$tddepresDist$NoBoot[iplot,]
          x <-   c(xx,rev(xx))
          y <- c(yy2,rev(yy1))
          x <- x[!is.na(y)   ]
          y <- y[!is.na(y)   ]
          polygon(x,y,col="grey88",border=NA)
          lines(xx,xx*0,col="grey20",lwd=2)
        }
        
        
    #    for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
    #    for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
      } 
    
      if ( ptype == 3){
        
        cyplot <- "CDFList[[isel]]$tddepresDist[istat,]"
        cxplot <- "CDFList[[isel]]$xseq"
        ylabel <- "Dew point depression [degC]"
        ylimits <- ylimforced 
        xlimits <- c(0,1)
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="cum. prob [0:1]")
        
        for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
        for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
      } 
      
      
      if ( ptype == 4){
        
        cyplot <- "CDFList[[isel]]$tdDist[iplot,]-CDFList[[1]]$tdDist[iplot,]"
        cxplot <- "CDFList[[isel]]$xseq"
        ylabel <- "Anom. dew point  wrt obs 1991-2020 [degC]"
        if (ytitle != "no")  { ylabel = cstat }
        ylimits <- ylimforced 
        xlimits <- c(0,1)
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="cum. prob [0:1]")
        
        grid()
        abline(h=0,col="grey80",lty=2)
        
        #    for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
        #    for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
        
        if (luncertainty){
          isel <- 1
          xx <- eval(parse(text=cxplot))
          yy2 <- CDFList[["uncertainty"]]$tdDist$P95[iplot,] - CDFList[["uncertainty"]]$tdDist$NoBoot[iplot,]
          yy1 <- CDFList[["uncertainty"]]$tdDist$P05[iplot,] - CDFList[["uncertainty"]]$tdDist$NoBoot[iplot,] 
          x <-   c(xx,rev(xx))
          y <- c(yy2,rev(yy1))
          x <- x[!is.na(y)   ]
          y <- y[!is.na(y)   ]
          polygon(x,y,col="grey88",border=NA)
          lines(xx,xx*0,col="grey20",lwd=2)
        }
        
      } 
    
      
      if ( ptype == 5){
        
        cyplot <- "CDFList[[isel]]$tdDist[istat,]"
        cxplot <- "CDFList[[isel]]$xseq"
        ylabel <- "Dew point [degC]"
        ylimits <- ylimforced 
        xlimits <- c(0,1)
        plot(NA,NA,xlim=xlimits,ylim=ylimits,ylab=ylabel,xlab="cum. prob [0:1]")
        
        for (hh in hlines) {abline(h=hh,col="grey80",lty=2)}
        for (vv in vlines) {abline(v=vv,col="grey80",lty=2)}
      } 
      
      pcol <- c()
      plegend <- c()
      
      if ( ptype <= 2 || ptype == 4) {
        for (isel in psel) {
          
          xx <- eval(parse(text=cxplot))
          yy <- eval(parse(text=cyplot))
          
          pcol <- c(pcol,EXPLIST[[isel]]$col)
          plegend <- c(plegend,EXPLIST[[isel]]$pname)
          
          lines(xx,yy,col=EXPLIST[[isel]]$col,lw=EXPLIST[[isel]]$lw)
        }
      }
     
      if ( ptype == 3 || ptype == 5 ) {
        csels <- c( "all","wet","90%", "99%", "99.9%" )
        ctmp <- c("black","green","cyan","blue", "magenta")
        
        for (csel in csels) {
          
          isel <- psel
          istat <- which(names(CDFList[[isel]]$prthres) == csel)
          
          xx <- eval(parse(text=cxplot))
          yy <- eval(parse(text=cyplot))
          pcol <- ctmp[csel == csels]
          ytmp <- format(CDFList[[isel]]$prthres[istat],digits = 3)
          plegend <- plegend <- c(plegend,paste0(csel," (",ytmp,")"))
           lines(xx,yy,col=pcol,lw=2)
        }
        
        pcol = ctmp
      }
      
      if ( ptext != "" ){
        if (ilegend == 0) { text(xlimits[2],ylimits[1],ptext )}
        if (ilegend == 1) { text(xlimits[1],ylimits[2],ptext,pos=4 )}
        if (ilegend == 2) { text(xlimits[2],ylimits[1],ptext,pos=2)}
      }
      
      if  ( ilegend != 0 ) {
        if (ilegend == 1) {legend(xlimits[2],ylimits[1],plegend,col=pcol,xjust=1,yjust=0,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
        if (ilegend == 2) {legend(xlimits[1],ylimits[2],plegend,col=pcol,xjust=0,yjust=1,cex=0.8,lw=2,seg.len=0.5,box.lwd=0) }
      }
    
      
      
      dev.off()
      
      
    }
    
    
    
    
    
    for (io in 1:2) {
      
      cdum <- def_olist[[io]]$cdum  
      psel <- def_olist[[io]]$psel
      cmodel <- names(def_olist)[io] 
      cdum <- paste0(cdum,ldash,seasonid)
      
      if (luncertainty){ 
          for (xvar in c("tddepresDist", "tdDist", "tasDist")){
           yy <- gl.UncertaintyBoot(ListIn=ResultsBootstrap$CDFlist,StatName = xvar)
            CDFDataList$uncertainty[[xvar]] <- yy[[xvar]]
          }
      }
      
      
      
      ptype = 1
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2])
      pptext = paste("50th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P50")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 2,ptext = pptext,ylimforced = c(0,9),pctl=0.5,luncertainty = luncertainty)
     
      pptext = paste("20th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P20")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 0,ptext = pptext,ylimforced = c(0,5),pctl=0.2,luncertainty = luncertainty)
      
      pptext = paste("80th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P80")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 0,ptext = pptext,ylimforced = c(0,12),pctl=0.8,luncertainty = luncertainty)
    
      pptext = paste("90th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P90")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 0,ptext = pptext,ylimforced = c(0,15),pctl=0.9, luncertainty = luncertainty)
    
     
      ptype = 2  # dew point
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2])
      pptext = paste("50th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P50_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 2,ptext = pptext,ylimforced = c(10,20),pctl=0.5, luncertainty = luncertainty)
      
      pptext = paste("20th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P20_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 0,ptext = pptext,ylimforced = c(6,16),pctl=0.2, luncertainty = luncertainty)
      
      pptext = paste("80th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P80_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 0,ptext = pptext,ylimforced = c(10,20),pctl=0.8, luncertainty = luncertainty)
      
      pptext = paste("90th percentile",ptext)
      cdum2 <- paste0(cdum,ldash,"CDFDPdep_P90_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdist(CDFDataList,EXPLIST,psel=psel,plotfile = plotfile,ptype=ptype,ilegend = 0,ptext = pptext,ylimforced = c(15,25),pctl=0.9, luncertainty = luncertainty)
    
      ptype = 2 
      cstat <- "all"
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2],cstat)
      
      cdum2 <- paste0(cdum,ldash,"CDFsingle_all_anom")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=psel[-1],
                               plotfile = plotfile,ptype=2,ilegend = 2,
                               ptext = ptext,ylimforced = c(-7,7),cstat = cstat, luncertainty = luncertainty)
      
      cstat <- "wet"
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2],cstat)
      cdum2 <- paste0(cdum,ldash,"CDFsingle_wet_anom")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=psel[-1],
                               plotfile = plotfile,ptype=2,ilegend = 2,
                               ptext = ptext,ylimforced = c(-7,7),cstat = cstat, luncertainty = luncertainty)
     
      cstat <- "99%"
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2],cstat)
      cdum2 <- paste0(cdum,ldash,"CDFsingle_P99_anom")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=psel[-1],
                               plotfile = plotfile,ptype=2,ilegend = 2,
                               ptext = ptext,ylimforced = c(-7,7),cstat = cstat, luncertainty = luncertainty)
      
      
      ptype = 4 
      cstat <- "all"
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2],cstat)
      
      cdum2 <- paste0(cdum,ldash,"CDFsingle_all_anom_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=psel[-1],
                               plotfile = plotfile,ptype=ptype,ilegend = 2,
                               ptext = ptext,ylimforced = c(-7,7),cstat = cstat, luncertainty = luncertainty)
      
      cstat <- "wet"
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2],cstat)
      cdum2 <- paste0(cdum,ldash,"CDFsingle_wet_anom_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=psel[-1],
                               plotfile = plotfile,ptype=ptype,ilegend = 2,
                               ptext = ptext,ylimforced = c(-7,7),cstat = cstat, luncertainty = luncertainty)
      
      cstat <- "99%"
      ptext <- paste(cmodel,strsplit(cobs2,"_")[[1]][2],cstat)
      cdum2 <- paste0(cdum,ldash,"CDFsingle_P99_anom_TD")
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=psel[-1],
                               plotfile = plotfile,ptype=ptype,ilegend = 2,
                               ptext = ptext,ylimforced = c(-7,7),cstat = cstat, luncertainty = luncertainty)
      
      
      }
    
    
    
    
    psels_single <- c(1:15)
    
    for (isel in psels_single){
      
      ptype = 3
      cdum2 <- paste0(names(EXPLIST)[isel],ldash,"CDFDPsingle_allstat",ldash,EXPLIST[[isel]]$prname,ldash,seasonid)
      ctmp <- paste(areaid,EXPLIST[[isel]]$pname)
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=isel,
                               plotfile = plotfile,ptype=ptype,ilegend = 2,
                               ptext = ctmp,ylimforced = c(0,20),cstat = cstat)
      
     
      ptype = 5
      cdum2 <- paste0(names(EXPLIST)[isel],ldash,"CDFDPsingle_allstat_TD",ldash,EXPLIST[[isel]]$prname,ldash,seasonid)
      ctmp <- paste(areaid,EXPLIST[[isel]]$pname)
      plotfile <- paste0(dirplots,cdum2,dash,areaid,"_ptype",ptype,".eps")
      gl.PlotCDFdistFullsingle(CDFDataList,EXPLIST,psel=isel,
                               plotfile = plotfile,ptype=ptype,ilegend = 2,
                               ptext = ctmp,ylimforced = c(5,25),cstat = cstat)
      
      
    }
    

}  # end over plotset
