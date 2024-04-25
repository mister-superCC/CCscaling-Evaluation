#
# converts bootstrapped objects to the old format without bootstrapping; note only without circulation clustering
# also does not use the low and lowX humidity classes but mid and low instead (needs to solved in postprocessing TBD)
# 20220804

hdir <- Sys.getenv("HOME")
#source(paste0(hdir,"/mysoft/R_lib/mylib.R"))
#source(paste0(hdir,"/mysoft/R_lib/ScalingLib.R"))
wdir <- paste0(hdir,"/analysis/PRINCIPLES/analysis_github/")
source(paste0(wdir,"/ScalingLibSmall.R"))

setwd(wdir)

library("ncdf4",         lib.loc=Rlibdir)
library("ncdf4.helpers", lib.loc=Rlibdir)
library("abind",         lib.loc=Rlibdir)
library("rlist",         lib.loc=Rlibdir)

plotsets = c(1,3,4)  # plotset 2 (all SFR data) is not processed and not needed anyway 
evalp <- "1991-2010"
pto <- 19
ptc <- 5
ptr <- 8
CPprname <- "pr_mean"
#CPprname <- "pr_sample"

dirdata <- "data_hourly_noboot/"
dirdata_boot <- "data_hourly_bootstrap/"
seasonid <- "MJJAS"

if ( ! dir.exists(dirdata)) { dir.create(dirdata)}


for (CPprname in c("pr_mean","pr_sample")){
for (plotset in plotsets){
      
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
      
      EXPLIST = list(
        
        # OBS, some room for different time periods, but not bootstrapped so three times copy of same data
        #OBS   =    list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "black",  pt=pto, lw = 3,  pname = "OBS 1991-2020"),
        #OBS2   =   list(period= "2011-2020", modelid = cobs,       prname = "precip", col = "grey40",   pt=pto, lw = 3, pname = "OBS 2011-2020"),
        #OBS3   =   list(period= "1991-2010", modelid = cobs,       prname = "precip", col = "grey80",   pt=pto, lw = 3,  pname = "OBS 1991-2010"),
        
        OBS   =    list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "black",  pt=pto, lw = 3,  pname = "OBS 1991-2020"),
        OBS2   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey40",   pt=pto, lw = 3, pname = "OBS 1991-2020"),
        OBS3   =   list(period= "1991-2020", modelid = cobs,       prname = "precip", col = "grey80",   pt=pto, lw = 3,  pname = "OBS 1991-2020"),
        
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
      
      
      # now try to convert to non bootstrapped data
      iexps = 1:15
      ccols <- NULL
      pchs <- NULL
        
      
      for (iexp in iexps){
        
        ccols <- c(ccols,EXPLIST[[iexp]]$col)
        pchs <-  c(pchs, EXPLIST[[iexp]]$pt)
        
        varnameP <- EXPLIST[[iexp]]$prname
        cobs <- EXPLIST[[iexp]]$modelid
        
        EnvBootstrap <- new.env()
        
        # bootstrap data format
        
        if (iexp %in% 1:3){ # obs data set
          
          ddid = paste(varnameP,cobs,seasonid,"1991-2020",sep="-")
       
        } else {
       
          ddid = paste(varnameP,EXPLIST[[iexp]]$modelid,areaid,seasonid,EXPLIST[[iexp]]$period,sep="-")
          
        }
        
      
        fileobj <- paste0(dirdata_boot,"/Cluster0_AllStatBootstrap_",ddid,".Rdata")
        
        
        load(file=fileobj, envir = EnvBootstrap)
        
        print(paste(fileobj,"exists", file.exists(fileobj)))
        
        
        
        # raw data format
        dataid <- "Cluster0_RHsens"
        Rdatafile <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
        
        #
        ScalingData_low = EnvBootstrap$ScalingList$raw$RHlow
        ScalingData_mid  = EnvBootstrap$ScalingList$raw$RHmid
        ScalingData_high = EnvBootstrap$ScalingList$raw$RHhigh
        
        save(file=Rdatafile,ScalingData_low,ScalingData_mid,ScalingData_high,seasonid,varnameP) 
        
        dataid <- "Cluster"
        Rdatafile2 <- paste0(dirdata,dataid,"_",ddid,"-alluseoneseas.Rdata")
        
        CDFlist <- list()
        CDFlist$all = EnvBootstrap$CDFlist$raw
        ExtremesList <- list()
        ExtremesList$all = EnvBootstrap$ExtremesList$raw
        PDFlist <- list()
        PDFlist$all = EnvBootstrap$PDFlist$raw
        ScalingList <- list()
        ScalingList$all = EnvBootstrap$ScalingList$raw
        ScalingListTas <- list()
        ScalingListTas$all = EnvBootstrap$ScalingListTas$raw
        
        
        save(file=Rdatafile2,CDFlist, ExtremesList, PDFlist, ScalingList, ScalingListTas, seasonid, varnameP)
        
        print(Rdatafile)
      #  print(Rdatafile2)
        
        rm(EnvBootstrap)
        
        # now convert bootstrap files back
        
      
      }
      
      
}
}









