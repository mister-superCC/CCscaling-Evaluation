#
# note T is used for dew point, Tas for temperature

hdir <- Sys.getenv("HOME")
#source(paste0(hdir,"/mysoft/R_lib/mylib.R"))
#source(paste0(hdir,"/mysoft/R_lib/ScalingLib.R"))
wdir <- paste0(hdir,"/analysis/PRINCIPLES/analysis_github/") # put here your working directory
# wdir <- paste0("/cnrm/mosca/USERS/cortese/Scripts/inpscript/in_Geert/METEOFRANCE_CODE_20231026/")
setwd(wdir)

source(paste0(wdir,"/ScalingLibSmall.R"))


library("ncdf4",         lib.loc=Rlibdir)
library("ncdf4.helpers", lib.loc=Rlibdir)
library("abind",         lib.loc=Rlibdir)
library("rlist",         lib.loc=Rlibdir)

#seasonid = "SHY"
#pthres <- 0.1 # 

pthress = c(0.1,-1) # absolute and conditional percentiles
pthress = c(-1)

nbootsizeFR = 10  # changed from 15 to 10 years to be consistent with other

for (seasonid in c("MJJAS")) {
# for (seasonid in c( "SHY", "JJA", "SON" ,"MJJASON" ,"DJF", "MAM")) {
  for (pthres in pthress){   
    
    # custum settings
    if (pthres == 0.1) {
        dirplots = "plots_hourly_bootstrap/"
        dirdata =  "data_hourly_bootstrap/"
        pquant <- c(0.9,0.99,0.999,0.9999)
    }
    
    if (pthres == -1) {
      dirplots = "plots_hourly_bootstrap_abs/"
      dirdata =  "data_hourly_bootstrap_abs/"
      pquant <- c(0.99,0.999,0.9999,0.99999)
    }
    
    if ( ! dir.exists(dirplots)) { dir.create(dirplots)}
    if ( ! dir.exists(dirdata))  { dir.create(dirdata)}
    
    ExcTres <- c(10,20,40,60) # array used for weight to fit to pquant
    bindef <- list("bin_width" = 2, "bin_min" = -10, "bin_max" = 30, "bin_step" = 1 )
    
    # some plotting things
    yrange1 <- c(0,120)
    yrange2 <- c(0.5,90)
    ylabel1 <- "hourly precipitation [mm/hour]"
    xrange1 <- c(1.,1e-7)
    
    
    PROJLIST <- list(
      OBS1 = list(PROJECT = "OBS", DOMAIN = "NL-STATIONS" ),          
      OBS2 = list(PROJECT = "OBS", DOMAIN = "MEDSEA-STATIONS" ) ,     
      OBS3 = list(PROJECT = "OBS", DOMAIN = "NoMEDSEA-STATIONS" ) ,   
      OBS4 = list(PROJECT = "OBS", DOMAIN = "SFR-STATIONS" )          
    )
    
    iproject_processed <- 1:length(PROJLIST)
    iproject_processed <- 2:4
    iproject_processed <- 1:3
    
    for (iproject in iproject_processed) {
      
      
      PROJECT <- PROJLIST[[iproject]]$PROJECT
      DOMAIN <- PROJLIST[[iproject]]$DOMAIN
      PROJDOM = paste(PROJECT,DOMAIN,sep="-")
      
      print(paste("processing PROJECT = ", PROJECT, DOMAIN))
     
      # some default settings 
      varnameT    <- "td2m"
      varnameTas  <- "t2m"
      pmul = 1.
      tadd = 0.
       
      # observations
      if(PROJDOM=="OBS-NL-STATIONS"){
        
    	  lOneDinput = TRUE
        maskfile = "no"
        maskid <- "NLstations"
    
        EXPLIST = list( Period1 = list(name="1991-2020",col="black", stam = "KNMI_20201124_hourly",
                                   maskdir = NA, varnameP  = "precip",
                                   ncindir =  c(paste0(hdir,"/analysis/OBSERVATIONS/NL-extended/")),
                                   yearsel = c(1991,2020),nbootsize=10) )
        
      } else if (PROJDOM=="OBS-SFR-STATIONS"){
        
        lOneDinput = TRUE
        maskfile <- "mask_SFR.nc"
        maskid <- "SFRstations"
        
        EXPLIST = list( Period1 = list(name="1991-2020",col="black", stam = "MF_data_1991-2020",
                        maskdir = NA,varnameP  = "precip",
                        ncindir =  c(paste0(hdir,"/LINKS/DataC/MF/")),
                        yearsel = c(1991,2020),nbootsize=nbootsizeFR) )
    
      } else if (PROJDOM=="OBS-MEDSEA-STATIONS") {
        
        lOneDinput = TRUE
        maskfile <- "mask_MedSea.nc"
        maskid <- "MEDSEAstations"
        
         EXPLIST = list( Period1 = list(name="1991-2020",col="black", stam = "MF_data_1991-2020",
                        maskdir = NA,varnameP  = "precip",
                        ncindir =  c(paste0(hdir,"/LINKS/DataC/MF/")),
                        yearsel = c(1991,2020),nbootsize=nbootsizeFR) )
                         
      } else if (PROJDOM=="OBS-NoMEDSEA-STATIONS") {
        
        lOneDinput = TRUE
        maskfile <- "mask_NoMedSea.nc"
        maskid <- "NoMEDSEAstations"
        
         EXPLIST = list( Period1 = list(name="1991-2020",col="black", stam = "MF_data_1991-2020",
                        maskdir = NA,varnameP  = "precip",
                        ncindir =  c(paste0(hdir,"/LINKS/DataC/MF/")),
                        yearsel = c(1991,2020),nbootsize=nbootsizeFR) )
    
      }
    
      
      
      
      # ------------------------------------------------------------------------------------------------------------------------------
      # start analysis
      # ------------------------------------------------------------------------------------------------------------------------------
     
      for (iperiod in 1:length(EXPLIST)) { 
    
        
        varnameP  = paste0(EXPLIST[[iperiod]]$varnameP)
     
        # read in model/obs data
        
        nquant <- length(pquant)
        
        ltest <- dir.exists(EXPLIST[[iperiod]]$ncindir)
        if (any(ltest)) {
          ncindir <- EXPLIST[[iperiod]]$ncindir[min(which(ltest))]
        } else {
          print("no valid directory for netcdf files found; ERROR !!")
        }  
        
        
        if ( ! maskfile == "no") { 
          if ( ! is.na(EXPLIST[[iperiod]]$maskdir )) {
            maskfile <- paste0(ncindir,EXPLIST[[iperiod]]$maskdir,"/mask_",maskid,".nc") 
          } else { 
            maskfile <- paste0(wdir,maskfile) 
          }
        }   
           
           
        if ( ! file.exists(maskfile) & ! maskfile == "no" ) {
          
            print(paste("maskfile does not exist",maskfile))
          
        } else {
      
             
            tas_infile <- paste0(ncindir,paste0(varnameTas,"_",EXPLIST[[iperiod]]$stam,"_3hr_",EXPLIST[[iperiod]]$name,".nc"))
            td_infile <-  paste0(ncindir,paste0(varnameT,"_",EXPLIST[[iperiod]]$stam,"_3hr_",EXPLIST[[iperiod]]$name,".nc"))
            pr_infile <-  paste0(ncindir,paste0(varnameP,"_",EXPLIST[[iperiod]]$stam,"_1hr_",EXPLIST[[iperiod]]$name,".nc"))
              
            if (PROJECT == "OBS") {
                tas_infile <- paste0(ncindir,paste0(EXPLIST[[iperiod]]$stam,".nc")) 
                td_infile  <- paste0(ncindir,paste0(EXPLIST[[iperiod]]$stam,".nc")) 
                pr_infile  <- paste0(ncindir,paste0(EXPLIST[[iperiod]]$stam,".nc"))
            }
            
          # read in data 
            indataP   <- gl.ReadFromNetcdf1D(ncinfile   = pr_infile, varnameP,   maskfile = maskfile)
            indataT   <- gl.ReadFromNetcdf1D(ncinfile   = td_infile, varnameT,   maskfile = maskfile)
            indataTas <- gl.ReadFromNetcdf1D(ncinfile  = tas_infile, varnameTas, maskfile = maskfile)
     
            timeP    <-   indataP$time
            timeT    <-   indataT$time
            timeTas  <- indataTas$time
            
            if (identical(timeT,timeTas)) {
              print("Time tas and td are identical; continue")
            } else {
              print("WARNING: Time tas and td are NOT identical; continue with care")
            }
        
            # first get season and year if present
        
            lsel <- month(timeP) %in% gl.monthlist(seasonid)     # mm <- as.numeric(format(timeP,"%m")) this is more robust one, works also with 360 day cal.
            if ( "yearsel" %in% names(EXPLIST[[iperiod]])) {
              yysel <- seq(EXPLIST[[iperiod]]$yearsel[1],EXPLIST[[iperiod]]$yearsel[2],1)
              lsel <- lsel & year(timeP) %in% yysel
            }
            pr_season = indataP$data[,lsel]
            timeP_season <- timeP[lsel]
          
            lsel <- month(timeT) %in% gl.monthlist(seasonid)
            if ( "yearsel" %in% names(EXPLIST[[iperiod]])) {
              yysel <- seq(EXPLIST[[iperiod]]$yearsel[1],EXPLIST[[iperiod]]$yearsel[2],1)
              lsel <- lsel & year(timeT) %in% yysel
            }
            td_season = indataT$data[,lsel]
            timeT_season <- timeT[lsel]
            tas_season <- indataTas$data[,lsel]
        
            # free up some memory
            rm(indataP,indataT,indataTas)
            
            # matching temp data to time precipitation # take hour minus 4 
            # iTimeTrans <- gl.TimeMatchClosest(timeP_season,timeT_season,6,4,12) # 4 hours before
            iTimeTrans <- gl.TimeMatchClosestN(timeP_season,timeT_season,4)
            td_season  <- td_season[,iTimeTrans]
            tas_season <- tas_season[,iTimeTrans]
            timeT_season <- timeT_season[iTimeTrans]
            
            # variables available pr_season, td_season, tas_season, timeP_season, timeT_season
            td_season = td_season + tadd
            tas_season = tas_season + tadd
             
            dataid = paste(varnameP,PROJDOM,maskid,seasonid,EXPLIST[[iperiod]]$name,sep="-")
    
            # ready for analysis
            # ----------------------------------------------------------------------------------------------------------------------------- 
            
            # create list to store the bootstrapped data
             
            ScalingList = list()
            ScalingListTas = list()
            ExtremesList = list()
            PDFlist = list()
            CDFlist <- list()
        
            # do here bootstrap !!!!!  to be done 
            nboot = 100
            ayears <- unique(year(timeP_season))
            nbsize = EXPLIST[[iperiod]]$nbootsize # sample size of bootstrap resample
                       # nbsize = length(ayears) # sample size of bootstrap resample
            
            for (iboot in 0:nboot) {
              
                  ncb <- paste0("iboot",iboot)
              
                  if ( nbsize/length(ayears) < 0.7) {
                    selyears <- sample(ayears,nbsize,replace=FALSE) 
                  } else {
                    selyears <- sample(ayears,nbsize,replace=TRUE) 
                  }
                    
                  if (iboot == 0) {
                    selyears = ayears
                    ncb <- "raw"
                    }
                  lsel <- year(timeP_season) %in% selyears
                  prsel <- pr_season[,lsel]
                  tdsel <- td_season[,lsel]
                  tassel <- tas_season[,lsel]
                  dpdepres <- tassel - tdsel
                  
                  ExtremesData <- gl.pr_extremes(prsel,nblock=24)
                  ScalingData     <- gl.Scaling(bindef=bindef,pquant = pquant, temp_paired = tdsel,  pr_paired = prsel, pthres=pthres )
                  ScalingDataTas  <- gl.Scaling(bindef=bindef,pquant = pquant, temp_paired = tassel, pr_paired = prsel, pthres=pthres )
                  PDFdata <- gl.PDF(prsel,ymax=200)
                
                  ScalingList[[ncb]]     <- ScalingData
                  ScalingListTas[[ncb]]  <- ScalingDataTas
                  ExtremesList[[ncb]]    <- ExtremesData
                  ExtremesList[[ncb]]$years <- selyears
                  PDFlist[[ncb]]         <- PDFdata
                  
                  # subselections based on surface dew point depresssion
                  lsel <- dpdepres < 3. #  high RH
                  ScalingData_high  <- gl.Scaling(bindef=bindef,pquant = pquant, temp_paired = tdsel[lsel], 
                                                  pr_paired = prsel[lsel], pthres=pthres )
                  ScalingList[[ncb]]$RHhigh  <- ScalingData_high
                  
                  lsel <- dpdepres >= 3. & dpdepres <= 6.# mid RH
                  ScalingData_mid  <- gl.Scaling(bindef=bindef,pquant = pquant, temp_paired = tdsel[lsel], 
                                                 pr_paired = prsel[lsel], pthres=pthres )
                  ScalingList[[ncb]]$RHmid  <- ScalingData_mid
                  
                  lsel <- dpdepres >= 6. # low RH
                  ScalingData_low  <- gl.Scaling(bindef=bindef,pquant = pquant, temp_paired = tdsel[lsel], 
                                                  pr_paired = prsel[lsel], pthres=pthres )
                  ScalingList[[ncb]]$RHlow  <- ScalingData_low
                  
                 
                  CDFobj <- gl.CDFexceedancePr(tdsel,tassel,prsel,narm=TRUE)
                  CDFlist[[ncb]] <- CDFobj
          
            
            } # bootstrap samples
        
            fileobj <- paste0(dirdata,"/Cluster0_AllStatBootstrap_",dataid,".Rdata")
            save(CDFlist,ScalingList,ScalingListTas,ExtremesList,PDFlist, seasonid, varnameP,version = 2,
                          file=fileobj )
               
               
        } # if mask file
            
        
     
      } # iperiod
      
    
#      plotfile <-  paste0(dirplots,"/Cluster0_AllStatBootstrap_",dataid,".eps") 
      plotfile <-  paste0(dirplots,"/Cluster0_AllStatBootstrap_",dataid) 
      gl.PlotScalingBoot(ScalingList,plotfile=plotfile)
        
    } # end loop of project
  

}  # end over seasonid and pthres
}





