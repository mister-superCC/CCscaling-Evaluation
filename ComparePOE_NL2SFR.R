#
# note T is used for dew point, Tas for temperature

hdir <- Sys.getenv("HOME")
wdir <- paste0(hdir,"/analysis/PRINCIPLES/analysis_github/")
source(paste0(wdir,"/ScalingLibSmall.R"))
iepsout = 1 # force eps file
setwd(wdir)

library("ncdf4",         lib.loc=Rlibdir)
library("ncdf4.helpers", lib.loc=Rlibdir)
library("abind",         lib.loc=Rlibdir)
library("rlist",         lib.loc=Rlibdir)

plotdir <- "plots_CDF_NLvsSFR/"
if (! dir.exists(plotdir)) dir.create(plotdir)

# custum settings

CPprname <- "pr_mean"
CPprname <- "pr_sample"
seasonid <- "MJJAS"

# areaidNL <-  "simpleNL.all.oroglt400" # "NLstations.land.all"
# cobsNL   <- "OBS-NL-STATIONS"
# areaidSFR <- "FRANCESOUTHbox.land.oroglt400" # "SFRstations.land.all"
# cobsSFR <- "OBS-SFR-STATIONS"
  
areaidNL <-  "NLstations.land.all"
cobsNL   <- "OBS-NL-STATIONS-NLstations"
dirdata_boot <-  "data_hourly_bootstrap/"

iopt = 2

if ( iopt == 1) {
    areaidSFR <-  "SFRstations.land.all"
    cobsSFR <- "OBS-SFR-STATIONS-SFRstations" 
    FRid <- "SFR"
}

if ( iopt == 2) {
    areaidSFR <-  "SFRstationsNoMedSea.land.all"
    cobsSFR <-   "OBS-NoMEDSEA-STATIONS-NoMEDSEAstations"
    FRid <- "SFR-cent"
}

if ( iopt == 3) {
    areaidSFR <-  "SFRstationsMedSea.land.all"
    cobsSFR <-   "OBS-MEDSEA-STATIONS-MEDSEAstations"
    FRid <- "SFR-med"
}




# directories in and output
dirdata =  "data_hourly_noboot/"


evalp <- "1991-2010"
pto <- 19
ptc <- 5
ptr <- 8


cobs <- cobsNL
EXPLIST_NL = list(
  
      # OBS 
      OBS_NL   =    list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "black",  pt=pto, lw = 3,  pname = "OBS 1991-2020"),
      OBS_NL2   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey40",   pt=pto, lw = 3, pname = "OBS 1991-2020"),
      OBS_NL1   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey80",   pt=pto, lw = 3,  pname = "OBS 1991-2020"),
      
      # CPMs
      ETH =      list(period= "CTL", modelid = "CPM-ERAINT-ETH-ALP3",         prname = CPprname, col = "red",    pt=ptc, lw = 1.5,  pname = "COSMO"),
      METO =     list(period= "CTL", modelid = "CPM-ERAINT-METO-ALP3",        prname = CPprname, col = "magenta",pt=ptc, lw =  1.5, pname = "UKMO-UM"),
      HCLIMALP3 =list(period= "CTL", modelid = "CPM-ERAINT-HCLIM-ALP3",       prname = CPprname, col = "blue",   pt=ptc, lw =  1.5, pname = "HCLIM-ALP"),
      HCLIMNWE3 =list(period= "CTL", modelid = "CPM-ERAINT-HCLIM-NWE3",       prname = CPprname, col = "cyan",   pt=ptc, lw =  1.5, pname = "HCLIM-NWE"),
      CNRM      =list(period= "CTL", modelid = "CPM-ERAINT-CNRM-NWE3",        prname = CPprname, col = "green",  pt=ptc, lw =  1.5, pname = "AROME41"),
      
      # CRCMs
      CLMCOM =   list(period= evalp, modelid = "RCM-ERAINT-CLMcom-ETH-COSMO", prname = "pr",     col = "red",   pt=ptr, lw = 2,    pname = "CLMcom"),
      HADRM =    list(period= evalp, modelid = "RCM-ERAINT-MOHC-HadREM3",     prname = "pr",     col = "magenta",   pt=ptr, lw = 2,pname = "HadRM3"),
      RACMO =    list(period= evalp, modelid = "RCM-ERAINT-KNMI-RACMO22E",    prname = "pr",     col = "blue",pt=ptr, lw = 2,      pname = "RACMO"),
      RCAO =     list(period= evalp, modelid = "RCM-ERAINT-SMHI-RCA4",        prname = "pr",     col = "cyan",   pt=ptr, lw = 2,   pname = "RCA4"),
      HIRHAM   = list(period= evalp, modelid = "RCM-ERAINT-DMI-HIRHAM5",      prname = "pr",     col = "green4" ,   pt=ptr, lw = 2,pname = "HIRHAM5"),
      REMO =     list(period= evalp, modelid = "RCM-ERAINT-GERICS-REMO2015",  prname = "pr",     col = "yellow2",   pt=ptr, lw = 2,pname = "REMO"),
      ALADIN =   list(period= evalp, modelid = "RCM-ERAINT-CNRM-ALADIN63",    prname = "pr",     col = "green",   pt=ptr, lw = 2,  pname = "ALADIN")
      
)



cobs <- cobsSFR
EXPLIST_SFR = list(
  
  # OBS 
  OBS_SFR   =    list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "black",    pt=pto, lw = 3,  pname = "OBS 1991-2020"),
  OBS_SFR2   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey40",   pt=pto, lw = 3,  pname = "OBS 1991-2020"),
  OBS_SFR1  =   list(period= "1991-2020", modelid = cobs,        prname = "precip", col = "grey80",   pt=pto, lw = 3,  pname = "OBS 1991-2020"),
  
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



nexps <- length(EXPLIST_NL)
iexps <-  1:nexps

ccols <-c()
pchs <- c()

ExtremesDataListNL <- list()
ScalingDataListNL <- list()
RHDataListNL <- list()
ExtremesDataListSFR <- list()
ScalingDataListSFR <- list()
RHDataListSFR <- list()


for (iexp in iexps){
  
# data from NL
  
    areaid <- areaidNL
  
    EXPLIST <- EXPLIST_NL  # things from EXPLIST
    expid <- names(EXPLIST)[iexp]
    ccols <- c(ccols,EXPLIST[[iexp]]$col)
    pchs <-  c(pchs, EXPLIST[[iexp]]$pt)
    varnameP <- EXPLIST[[iexp]]$prname
    
    tmpenv <- new.env()
    

    if (iexp %in% 1:3){ # obs data set
      
      ddid = paste(varnameP,cobsNL,seasonid,"1991-2020",sep="-")
      
    } else {
      
      ddid = paste(varnameP,EXPLIST[[iexp]]$modelid,areaid,seasonid,EXPLIST[[iexp]]$period,sep="-")
      
    }
    
    dataid <- "Cluster0_RHsens"
    Rdatafile <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
    
    dataid <- "Cluster"
    Rdatafile2 <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
    
    
    
    load(Rdatafile,envir=tmpenv)
    RHDataListNL[[expid]] <- list("RHlow"  = tmpenv$ScalingData_lowX, 
                                  "RHmid"  = tmpenv$ScalingData_low,
                                  "RHhigh" = tmpenv$ScalingData_high )
    rm(tmpenv)
    
    tmpenv <- new.env()
    load(Rdatafile2,envir=tmpenv)
    ExtremesDataListNL[[expid]]      <- tmpenv$ExtremesList$all
    ExtremesDataListNL[[expid]]$CDF  <- tmpenv$CDFlist
    
    ScalingDataListNL[[expid]]  <- tmpenv$ScalingList$all 
    
    rm(tmpenv)
    
    
#    Data from SFR
    
    
    areaid <- areaidSFR

    EXPLIST <- EXPLIST_SFR  # things from EXPLIST
    expid <- names(EXPLIST)[iexp]
#    ccols <- c(ccols,EXPLIST[[iexp]]$col)
#    pchs <-  c(pchs, EXPLIST[[iexp]]$pt)
    varnameP <- EXPLIST[[iexp]]$prname
    
    
    if (iexp %in% 1:3){ # obs data set
      
      ddid = paste(varnameP,cobsSFR,seasonid,"1991-2020",sep="-")
      
    } else {
      
      ddid = paste(varnameP,EXPLIST[[iexp]]$modelid,areaid,seasonid,EXPLIST[[iexp]]$period,sep="-")
      
    }
    
    
    dataid <- "Cluster0_RHsens"
    Rdatafile <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
    
    dataid <- "Cluster"
    Rdatafile2 <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
    
    
    tmpenv <- new.env()
    load(Rdatafile,envir=tmpenv)
    
    RHDataListSFR[[expid]] <- list("RHlow"  = tmpenv$ScalingData_lowX, 
                                   "RHmid"  = tmpenv$ScalingData_low,
                                   "RHhigh" = tmpenv$ScalingData_high )
    rm(tmpenv)
    tmpenv <- new.env()
    load(Rdatafile2,envir=tmpenv)
    ExtremesDataListSFR[[expid]]      <- tmpenv$ExtremesList$all
    ExtremesDataListSFR[[expid]]$CDF  <- tmpenv$CDFlist
    
    ScalingDataListSFR[[expid]]  <- tmpenv$ScalingList$all 
    
    rm(tmpenv)
    
    
}

def_olist <- list(
  "CPMs" = list("psel" = c(1,4:8),
                "cdum" =  paste0("CPM_",CPprname)),
  "RCMs" = list("psel" = c(1,9:15),
                "cdum" =  "RCM_pr"))





ResultsBootstrapNL <- new.env()

varnameP <- "precip"  # OBS-NL-STATIONS-NLstations-SHY-1991-2020.Rdata*
dataid = paste(varnameP,"OBS-NL-STATIONS-NLstations",seasonid,"1991-2020",sep="-")
fileobj <- paste0(dirdata_boot,"/Cluster0_AllStatBootstrap_",dataid,".Rdata")
load(file=fileobj, envir = ResultsBootstrapNL)
 

ResultsBootstrapSFR <- new.env()

varnameP <- "precip"  # OBS-NL-STATIONS-NLstations-SHY-1991-2020.Rdata*
dataid = paste(varnameP,cobsSFR,seasonid,"1991-2020",sep="-")
fileobj <- paste0(dirdata_boot,"/Cluster0_AllStatBootstrap_",dataid,".Rdata")
load(file=fileobj, envir = ResultsBootstrapSFR)


#===================================================================================================
gl.UncertaintyBoot <- function(ListIn,StatName,nboot=100,LevelName = NA, LevelName2 = NA, rtype = "rel") {
  
  # ListIn <- ScalingListObsFit
  # StatName <- "xLOESS"
  
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

ExtListNL <- ResultsBootstrapNL$ExtremesList
ExtListSFR <- ResultsBootstrapSFR$ExtremesList

xx    <- ExtListNL$raw$pctls_extremes
yyNL  <- gl.UncertaintyBoot(ListIn=ExtListNL,StatName = "pvar_extremes")
yySFR <- gl.UncertaintyBoot(ListIn=ExtListSFR,StatName = "pvar_extremes")

RespList <- list()
for (cname in names(ExtListNL)){
  RespList[[cname]]$NL <- ExtListNL[[cname]]
  RespList[[cname]]$SFR <- ExtListSFR[[cname]]
}

RespBootstrap <- gl.UncertaintyBoot(ListIn=RespList,StatName = "pvar_extremes",LevelName = "SFR", LevelName2 = "NL")  

#

gl.PlotScalingDuo(ScalingDataListNL$OBS_NL,ScalingDataListSFR$OBS_SFR)
gl.PlotScalingDuo(RHDataListNL$OBS_NL$RHmid,RHDataListSFR$OBS_SFR$RHmid)
gl.PlotScalingDuo(RHDataListNL$OBS_NL$RHhigh,RHDataListSFR$OBS_SFR$RHhigh)
#gl.PlotScalingDuo(RHDataListNL$OBS_NL$RHlow,RHDataListSFR$OBS_SFR$RHlow)


FRid2 <- strsplit(cobsSFR,"-")[[1]][2]

fileout <- paste0(plotdir,def_olist[[1]]$cdum,"_",FRid2,"vsNL_PoE",ldash,seasonid,".eps")
source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
options(scipen=2)
par(cex=1.4)

ii=1
pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes / ExtremesDataListNL[[ii]]$pvar_extremes
xx <- 1-ExtremesDataListNL[[ii]]$pctls_extremes

pdiff = 100*(pdiff-1.)

pdiff_boot <- 100.*(ExtremesDataListSFR[[ii]]$pvar_extremes_boot / ExtremesDataListNL[[ii]]$pvar_extremes_boot - 1)
pdiff_P5_P95 <- apply(pdiff_boot,MARGIN=1,FUN=quantile,probs=c(0.05,0.5,0.95),na.rm=TRUE)

pdiff_P5_P95[1,] <- RespBootstrap$pvar_extremes$P05
pdiff_P5_P95[2,] <- RespBootstrap$pvar_extremes$P50
pdiff_P5_P95[3,] <- RespBootstrap$pvar_extremes$P95


plot(xx,pdiff+NA,col=ccols[1],type="l",lwd = 2,log="x",ylim=c(-50,100),xlim=c(0.1,5e-7),
     xlab = "probability of exceedance",ylab=paste("difference rainfall",FRid,"compared to NL [%]"))
grid()

# lines(xx,pdiff_P5_P95[2,],col="brown",lwd=2)
# lines(xx,pdiff_P5_P95[1,],col="brown",lwd=2)
# lines(xx,pdiff_P5_P95[3,],col="brown",lwd=2)

x <-   c(xx,rev(xx))
y <- c(pdiff_P5_P95[1,],rev(pdiff_P5_P95[3,]))

#y [x>0.5 | x < 10/(ExtremesData1$nall) | x2 < 0.5 ] <- NA
#             y [x>0.1 | x < 10/(ExtremesData1$nall) ] <- NA

y [x>0.04] <- NA
x <- x[!is.na(y)   ]
y <- y[!is.na(y)   ]
polygon(x,y,col="grey90",border=NA)


for (ii in def_olist$CPMs$psel) {
  pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes / ExtremesDataListNL[[ii]]$pvar_extremes
  pdiff = 100*(pdiff-1.)
  lines(xx,pdiff,col=ccols[ii],lwd = 2) }

text(1e-1,95,"CPMs",pos=4)
dev.off()







fileout <- paste0(plotdir,def_olist[[2]]$cdum,"_",FRid2,"vsNL_PoE",ldash,seasonid,".eps")
source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
options(scipen=2)
par(cex=1.4)

ii=1
pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes / ExtremesDataListNL[[ii]]$pvar_extremes
pdiff = 100*(pdiff-1.)
xx <- 1-ExtremesDataListNL[[ii]]$pctls_extremes


pdiff_boot <- 100.*(ExtremesDataListSFR[[ii]]$pvar_extremes_boot / ExtremesDataListNL[[ii]]$pvar_extremes_boot - 1)
pdiff_P5_P95 <- apply(pdiff_boot,MARGIN=1,FUN=quantile,probs=c(0.05,0.5,0.95),na.rm=TRUE)

pdiff_P5_P95[1,] <- RespBootstrap$pvar_extremes$P05
pdiff_P5_P95[2,] <- RespBootstrap$pvar_extremes$P50
pdiff_P5_P95[3,] <- RespBootstrap$pvar_extremes$P95


plot(xx,pdiff,col=ccols[1],type="l",lwd = 2,log="x",ylim=c(-50,100),xlim=c(0.1,5e-7),
     xlab = "probability of exceedance",ylab=paste("difference rainfall",FRid,"compared to NL [%]"))
grid()

#lines(xx,pdiff_P5_P95[2,],col="brown",lwd=2)
#lines(xx,pdiff_P5_P95[1,],col="brown",lwd=2)
#lines(xx,pdiff_P5_P95[3,],col="brown",lwd=2)

x <-   c(xx,rev(xx))
y <- c(pdiff_P5_P95[1,],rev(pdiff_P5_P95[3,]))

#y [x>0.5 | x < 10/(ExtremesData1$nall) | x2 < 0.5 ] <- NA
#             y [x>0.1 | x < 10/(ExtremesData1$nall) ] <- NA

y [x>0.04] <- NA
x <- x[!is.na(y)   ]
y <- y[!is.na(y)   ]
polygon(x,y,col="grey90",border=NA)


for (ii in def_olist$RCMs$psel) {
  pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes / ExtremesDataListNL[[ii]]$pvar_extremes
  pdiff = 100*(pdiff-1.)
  lines(xx,pdiff,col=ccols[ii],lwd = 2) }

text(1e-1,95,"RCMs",pos=4)
dev.off()



fileout <- paste0(plotdir,def_olist[[1]]$cdum,"_",FRid2,"vsNL_PoE_cond",ldash,seasonid,".eps")
source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
options(scipen=2)
par(cex=1.4)

ii=1
pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes_cond / ExtremesDataListNL[[ii]]$pvar_extremes_cond
xx <- 1-ExtremesDataListNL[[ii]]$pctls_extremes

plot(xx,pdiff,col=ccols[1],type="l",lwd = 2,log="x",ylim=c(0.5,2.5),xlim=c(0.1,5e-7),
     xlab = "probability of exceedance",ylab="fractional difference SFR/NL")
grid()

for (ii in def_olist$CPMs$psel) {
  pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes_cond / ExtremesDataListNL[[ii]]$pvar_extremes_cond
  lines(xx,pdiff,col=ccols[ii],lwd = 2) }

text(1e-1,1.95,"CPMs",pos=4)
dev.off()



fileout <- paste0(plotdir,def_olist[[2]]$cdum,"_",FRid2,"vsNL_PoE_cond",ldash,seasonid,".eps")
source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
options(scipen=2)
par(cex=1.4)

ii=1
pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes_cond / ExtremesDataListNL[[ii]]$pvar_extremes_cond
xx <- 1-ExtremesDataListNL[[ii]]$pctls_extremes

plot(xx,pdiff,col=ccols[1],type="l",lwd = 2,log="x",ylim=c(0.5,2.5),xlim=c(0.1,5e-7),
     xlab = "probability of exceedance",ylab="fractional difference SFR/NL")
grid()

for (ii in def_olist$RCMs$psel) {
  pdiff <- ExtremesDataListSFR[[ii]]$pvar_extremes_cond / ExtremesDataListNL[[ii]]$pvar_extremes_cond
  lines(xx,pdiff,col=ccols[ii],lwd = 2) }

text(1e-1,1.95,"RCMs",pos=4)
dev.off()



