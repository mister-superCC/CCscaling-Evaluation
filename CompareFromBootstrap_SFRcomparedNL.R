#
# note T is used for dew point, Tas for temperature

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
# custum settings

CPprname <- "pr_mean"
CPprname <- "pr_sample"
seasonid <- "MJJAS"
dirdata <-  "data_hourly_bootstrap/"
plotdir <- "plots_comparisonSFR2NL/"

cptype <- "abs"
cptype = "rel"

if (cptype == "abs"){
  dirdata <-  "data_hourly_bootstrap_abs/"
  plotdir <- "plots_comparisonSFR2NL_abs/"
}

if ( ! dir.exists(plotdir)) dir.create(plotdir)

areaid1 <-  "NLstations.land.all"
cobs1 <-   "NL-STATIONS-NLstations"
cobsid1 <- "NL"

iopt = 2
if ( iopt == 1) {
  areaid2 <-  "SFRstations.land.all"
  cobs2 <- "SFR-STATIONS-SFRstations" 
  cobsid2 <- "SFR"
}
if ( iopt == 2) {
  areaid2 <-  "SFRstationsNoMedSea.land.all"
  cobs2 <-   "NoMEDSEA-STATIONS-NoMEDSEAstations"
  cobsid2 <- "SFR-cent"
}
if ( iopt == 3) {
  areaid2 <-  "SFRstationsMedSea.land.all"
  cobs2 <-   "MEDSEA-STATIONS-MEDSEAstations"
  cobsid2 <- "SFR-med"
}


evalp <- "1991-2010"
pto <- 19
ptc <- 5
ptr <- 8


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



EXPLIST = list(
  
  # OBS 
  OBS   =    list(period= "1991-2020", modelid = "OBS",       prname = "precip", col = "black",  pt=pto, lw = 3,  pname = "OBS 1991-2020"),
  
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



nexps <- length(EXPLIST)
iexps <-  1:nexps
#iexps <-  1:1

ccols <-c()
pchs <- c()

AllStat1 <-  list()
AllStat2 <-  list()


for (iexp in iexps){
  
  #  Cluster0_AllStatBootstrap_precip-OBS-SFR-STATIONS-SFRstations-SHY-1991-2020.Rdata
  ResultsBootstrap1 <- new.env()
  
  varnameP <- EXPLIST[[iexp]]$prname 
  modelid <- EXPLIST[[iexp]]$modelid
  periodid <- EXPLIST[[iexp]]$period
  expname <- names(EXPLIST)[iexp]
  
  if (names(EXPLIST)[[iexp]] == "OBS") {
    dataid = paste(varnameP,modelid,cobs1,seasonid,periodid,sep="-")
  } else {
    dataid = paste(varnameP,modelid,areaid1,seasonid,periodid,sep="-")
  }
  fileobj <- paste0(dirdata,"Cluster0_AllStatBootstrap_",dataid,".Rdata")
  load(file=fileobj, envir = ResultsBootstrap1)
  AllStat1[[expname]] <- ResultsBootstrap1
  
  ResultsBootstrap2 <- new.env()
  
  if (names(EXPLIST)[[iexp]] == "OBS") {
    dataid = paste(varnameP,modelid,cobs2,seasonid,periodid,sep="-")
  } else {
    dataid = paste(varnameP,modelid,areaid2,seasonid,periodid,sep="-")
  }
  fileobj <- paste0(dirdata,"Cluster0_AllStatBootstrap_",dataid,".Rdata")
  load(file=fileobj, envir = ResultsBootstrap2)
  AllStat2[[expname]] <- ResultsBootstrap2
  
}

# AllStat1 and AllStat2 contain all postprocessing data ; a lot !

for (iexp in iexps){
  
  modelid <- EXPLIST[[iexp]]$modelid
  varnameP <- EXPLIST[[iexp]]$prname 
  periodid <- EXPLIST[[iexp]]$period
  modelname <- EXPLIST[[iexp]]$pname
  
  dataid = paste(varnameP,modelid,cobsid1,cobsid2,seasonid,periodid,sep="-")
  plotfile = paste0(plotdir,"/ScalingTD",ldash,dataid,".eps")
  gl.PlotScalingDuo(AllStat1[[iexp]]$ScalingList$raw,AllStat2[[iexp]]$ScalingList$raw,plotfile = plotfile,
                    plottitlelt = c(paste(modelname,cobsid1),paste(modelname,cobsid2)))
  
  
  Scaling1 <- list()
  Scaling2 <- list()
  plotfile = paste0(plotdir,"/ScalingTD",ldash,dataid)
  
  for (ctmp in names(AllStat1[[1]]$ScalingList)) {
    Scaling1[[ctmp]] <- AllStat1[[iexp]]$ScalingList[[ctmp]]
    Scaling2[[ctmp]] <- AllStat2[[iexp]]$ScalingList[[ctmp]]
  }
  gl.PlotScalingDuoBoot(Scaling1,Scaling2,plotfile = plotfile,
                        plottitlelt = c(paste(modelname,cobsid1),paste(modelname,cobsid2)))
  

  
  plotfile = paste0(plotdir,"/ScalingTD_RHhigh",ldash,dataid,".eps")
  gl.PlotScalingDuo(AllStat1[[iexp]]$ScalingList$raw$RHhigh,AllStat2[[iexp]]$ScalingList$raw$RHhigh,plotfile = plotfile,
                    plottitlelt = c(paste(modelname,cobsid1),paste(modelname,cobsid2)))

  
  Scaling1 <- list()
  Scaling2 <- list()
  plotfile = paste0(plotdir,"/ScalingTD_RHhigh",ldash,dataid)
  
  for (ctmp in names(AllStat1[[1]]$ScalingList)) {
    Scaling1[[ctmp]] <- AllStat1[[iexp]]$ScalingList[[ctmp]]$RHhigh
    Scaling2[[ctmp]] <- AllStat2[[iexp]]$ScalingList[[ctmp]]$RHhigh
  }
  gl.PlotScalingDuoBoot(Scaling1,Scaling2,plotfile = plotfile,
                        plottitlelt = c(paste(modelname,cobsid1),paste(modelname,cobsid2)))
  
  
  
  plotfile = paste0(plotdir,"/ScalingTD_RHlow",ldash,dataid,".eps")
  gl.PlotScalingDuo(AllStat1[[iexp]]$ScalingList$raw$RHlow,AllStat2[[iexp]]$ScalingList$raw$RHlow,plotfile = plotfile,
                    plottitlelt = c(paste(modelname,cobsid1),paste(modelname,cobsid2)))
  
  
  Scaling1 <- list()
  Scaling2 <- list()
  plotfile = paste0(plotdir,"/ScalingTD_RHlow",ldash,dataid)
  
  for (ctmp in names(AllStat1[[1]]$ScalingList)) {
    Scaling1[[ctmp]] <- AllStat1[[iexp]]$ScalingList[[ctmp]]$RHlow
    Scaling2[[ctmp]] <- AllStat2[[iexp]]$ScalingList[[ctmp]]$RHlow
  }
  gl.PlotScalingDuoBoot(Scaling1,Scaling2,plotfile = plotfile,
                        plottitlelt = c(paste(modelname,cobsid1),paste(modelname,cobsid2)))
  
  
}











#source(paste0(hdir,"/mysoft/R_lib/ScalingLib.R"))
#gl.PlotScalingBoot(AllStat$HIRHAM$ScalingList,linterpolate = TRUE,rlspan=0.5)


ListIn1 <- AllStat1$ OBS$CDFlist
ListIn2 <- AllStat2$ OBS$CDFlist
varname <- "tddepresDist"
#varname <- "tdDist"
xx  <- gl.UncertaintyBoot(ListIn=ListIn1,StatName = "xseq")$xseq$NoBoot
yy1 <- gl.UncertaintyBoot(ListIn=ListIn1,StatName = varname)[[varname]]$NoBoot
yy2 <- gl.UncertaintyBoot(ListIn=ListIn2,StatName = varname)[[varname]]$NoBoot


ListIn1 <- AllStat1$OBS$ScalingList
ListIn2 <- AllStat2$OBS$ScalingList
ctmps <- names(ListIn1)
xx1 <- AllStat1$OBS$ScalingList$raw$xvar_mean
xx2 <- AllStat2$OBS$ScalingList$raw$xvar_mean


#===================================================================================================
gl.PlotLocal <- function(xx1,xx2,yy1,yy2,xlab="dew point temperature",ylab="frequency",
                         xlim=c(4,24),ylim=range(yy1$NoBoot,na.rm=TRUE),plotfile = "test") {
  
  iepsout <- 5
  fileout <- plotfile
  source(paste0(hdir,"/mysoft/R_lib/initeps.R"),local=TRUE)
  options(scipen=2)
  par(cex=1.4)
  
  plot(xx1,yy1$NoBoot+NA,type="l",lwd=2,xlab=xlab,ylab = ylab,xlim=xlim,ylim=ylim)
  x <-   c(xx1,rev(xx1))
  y <- c(yy1$P05,rev(yy1$P95))
  x <- x[!is.na(y)   ]
  y <- y[!is.na(y)   ]
  polygon(x,y,col=rgb(0,0,0,0.1),border=NA)
  lines(xx1,yy1$NoBoot,type="l",lwd=2,col="black")
  x <-   c(xx2,rev(xx2))
  y <- c(yy2$P05,rev(yy2$P95))
  x <- x[!is.na(y)   ]
  y <- y[!is.na(y)   ]
  polygon(x,y,col=rgb(1,0,0,0.1),border=NA)
  lines(xx2,yy2$NoBoot,type="l",lwd=2,col="red")
  legend(xlim[1]+1,0.95*ylim[2],c(cobsid1,cobsid2),col=c("black","red"),xjust=0,yjust=1,cex=0.8,lw=2,seg.len=0.5,box.lwd=0)
  dev.off()
}



# wet hour frequency; all events

TmpList1 <- list()
TmpList2 <- list()
for (ctmp in ctmps){
  ntmp <- ListIn1[[ctmp]]$Nall
  tmp <- ListIn1[[ctmp]]$Nwet/ListIn1[[ctmp]]$Nall
  tmp[ntmp<100] <- NA
  TmpList1[[ctmp]]$fwet  <- tmp
  
  ntmp <- ListIn2[[ctmp]]$Nall
  tmp <- ListIn2[[ctmp]]$Nwet/ListIn2[[ctmp]]$Nall
  tmp[ntmp<100] <- NA
  TmpList2[[ctmp]]$fwet  <- tmp
} 
varname = "fwet"
yy1 <- gl.UncertaintyBoot(ListIn=TmpList1,StatName = varname)[[varname]]
yy2 <- gl.UncertaintyBoot(ListIn=TmpList2,StatName = varname)[[varname]]


dataid <- paste(cobsid1,cobsid2,seasonid,sep="-")
gl.PlotLocal(xx1,xx2,yy1,yy2,plotfile=paste0(plotdir,"FreqAll_",dataid),ylim=c(0,0.2))
# ploting



csubsel <- "RHmid"
TmpList1 <- list()
TmpList2 <- list()
for (ctmp in ctmps){
  ntmp <- ListIn1[[ctmp]][[csubsel]]$Nall
  tmp <- ListIn1[[ctmp]][[csubsel]]$Nwet/ListIn1[[ctmp]][[csubsel]]$Nall
  tmp[ntmp<100] <- NA
  TmpList1[[ctmp]]$fwet  <- tmp
  
  ntmp <- ListIn2[[ctmp]][[csubsel]]$Nall
  tmp <- ListIn2[[ctmp]][[csubsel]]$Nwet/ListIn2[[ctmp]][[csubsel]]$Nall
  tmp[ntmp<100] <- NA
  TmpList2[[ctmp]]$fwet  <- tmp
} 
varname = "fwet"
yy1 <- gl.UncertaintyBoot(ListIn=TmpList1,StatName = varname)[[varname]]
yy2 <- gl.UncertaintyBoot(ListIn=TmpList2,StatName = varname)[[varname]]

dataid <- paste(cobsid1,cobsid2,seasonid,sep="-")
gl.PlotLocal(xx1,xx2,yy1,yy2,plotfile=paste0(plotdir,"Freq_",csubsel,ldash,dataid),ylim=c(0,0.2))


