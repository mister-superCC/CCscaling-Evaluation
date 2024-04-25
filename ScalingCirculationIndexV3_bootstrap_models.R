#
# note T is used for dew point, Tas for temperature

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

lRdatreprocess <-  FALSE

pthres <- -1  # pthres <- 0.1 #
#pthres <- 0.1  # pthres <- 0.1 #

# custum settings
dirplots = "plots_hourly_bootstrap/"
dirdata =  "data_hourly_bootstrap/"
pquant <- c(0.9,0.99,0.999,0.9999)

# link to the input data containg the R object files
rodir =  paste0(hdir,"/dropbox/tmp_data/DATA_SCALING_PRINCIPLES/")


if (pthres == -1) {
  pquant <- c(0.99,0.999,0.9999,0.99999)
  dirplots = "plots_hourly_bootstrap_abs/"
  dirdata =  "data_hourly_bootstrap_abs/"
  }

rcmindir <- "/nobackup_1/users/lenderin/EURO-CORDEX/PRINCIPLES/resized_10-year-files/"

if ( ! dir.exists(dirplots)) { dir.create(dirplots)}
if ( ! dir.exists(dirdata))  { dir.create(dirdata)}

ExcTres <- c(10,20,40,60) # array used for weight to fit to pquant
bindef <- list("bin_width" = 2, "bin_min" = -10, "bin_max" = 30, "bin_step" = 1 )

# some plotting things
yrange1 <- c(0,120)
yrange2 <- c(0.5,90)
ylabel1 <- "hourly precipitation [mm/hour]"
xrange1 <- c(1.,1e-7)

maskids   <-  c( "SFRstationsNoMedSea.land.all" , "SFRstationsMedSea.land.all", "NLstations.land.all" )
# seasonids <-  c( "SHY", "JJA", "SON" ,"MJJASON" ,"DJF", "MAM")
seasonids <-  c( "MJJAS" )


for (maskid in maskids) {
  for (seasonid in seasonids) {


      # reads input if run in script
      if (exists("loadinfile")) {
        if (file.exists(loadinfile)) {
          source(loadinfile)
        } else {
          print(paste(loadinfile," DOES NOT EXIST & NOT READ IN"))
        }
      }


      PROJLIST <- list(
        RC1 = list(PROJECT = "RCM", DOMAIN="ERAINT-DMI-HIRHAM5"),     # principles runs 1
        RC2 = list(PROJECT = "RCM", DOMAIN="ERAINT-KNMI-RACMO22E" ), # 2
        RC3 = list(PROJECT = "RCM", DOMAIN="ERAINT-SMHI-RCA4" ),
        RC4 = list(PROJECT = "RCM", DOMAIN="ERAINT-GERICS-REMO2015" ), # 4
        RC5 = list(PROJECT = "RCM", DOMAIN="ERAINT-MOHC-HadREM3" ),
        RC6 = list(PROJECT = "RCM", DOMAIN="ERAINT-CLMcom-ETH-COSMO"), #6
        RC7 = list(PROJECT = "RCM", DOMAIN="ERAINT-CNRM-ALADIN63"),
        CP1 = list(PROJECT = "CPM", DOMAIN="ERAINT-HCLIM-ALP3" ),     # 8 hclim runs
        CP2 = list(PROJECT = "CPM", DOMAIN="ERAINT-HCLIM-NWE3" ),
        CPETH  = list(PROJECT = "CPM", DOMAIN="ERAINT-ETH-ALP3" ),   #10
        CPMETO = list(PROJECT = "CPM", DOMAIN="ERAINT-METO-ALP3" ),
        CPCNRM = list(PROJECT = "CPM", DOMAIN="ERAINT-CNRM-NWE3" ) # 12
      )

      iproject_processed <- 1:length(PROJLIST)
#      iproject_processed <- 1

      for (iproject in iproject_processed) {


        PROJECT <- PROJLIST[[iproject]]$PROJECT
        DOMAIN <- PROJLIST[[iproject]]$DOMAIN
        PROJDOM = paste(PROJECT,DOMAIN,sep="-")

        print(paste("processing PROJECT = ", PROJECT, DOMAIN))

        # some default settings
        varnameT  <- "tdps"
        varnameTas  <- "tas"
        pmul = 3600.        # multiplication factor precip.
        tadd = -273.15      # addition to temperature

        # observations
        if(PROJDOM=="RCM-ERAINT-DMI-HIRHAM5"){
          EXPLIST = list(
              P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v1",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"DMI-HIRHAM5/NLcatted/")))
          )
        }

        if(PROJDOM=="RCM-ERAINT-KNMI-RACMO22E"){
          EXPLIST = list(
            P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22E_v1",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"KNMI-RACMO22E/NLcatted/")))
          )
        }

        if(PROJDOM=="RCM-ERAINT-SMHI-RCA4"){
          EXPLIST = list(
            P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"SMHI-RCA4/NLcatted/")))
          )
        }

        if(PROJDOM=="RCM-ERAINT-GERICS-REMO2015"){
          EXPLIST = list(
            P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v2",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"GERICS-REMO2015/NLcatted/")))
          )
          pmul = 1 # multiplication factor precip.

        }


        if(PROJDOM=="RCM-ERAINT-MOHC-HadREM3"){
          EXPLIST = list(
            P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_MOHC-HadREM3-GA7-05_v1",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"MOHC-HadREM3-GA7-05/NLcatted/")))
          )

        }

        if(PROJDOM=="RCM-ERAINT-CLMcom-ETH-COSMO"){
          EXPLIST = list(
             P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"CLMcom-ETH-COSMO-crCLIM-v1-1/NLcatted/")))
          )

        }

        # CNRM-ALADIN63

        if(PROJDOM=="RCM-ERAINT-CNRM-ALADIN63"){
          EXPLIST = list(
            P1991_2010 = list(name="1991-2010",col="black", stam = "EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CNRM-ALADIN63_v1",
                              maskdir = "masks/", varnameP  = "pr",
                              ncindir =  c(paste0(rcmindir,"CNRM-ALADIN63/NLcatted/")))
          )

        }

        # CPM runs

        if(PROJDOM=="CPM-ECEARTH-HCLIM-ALP3"){
          EXPLIST = list(
            CTL=list(name="CTL",col="black", stam = "his.CXFPS025.HCLIM38h1_CXFPS_fRACMOfECEARTH_r14_hist.AYEAR.1hr.G5x5",
                     maskdir = paste0(hdir,"/analysis/CXFPS/PREPROCESS/NC/"), varnameP  = "pr_mean",
                     ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r14_hist_ReGrid/"),
                                  paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r14_hist_ReGrid/")))
            ,
            CTL2 =list(name="CTL",col="black", stam = "his.CXFPS025.HCLIM38h1_CXFPS_fRACMOfECEARTH_r14_hist.AYEAR.1hr.G5x5",
                       maskdir = paste0(hdir,"/analysis/CXFPS/PREPROCESS/NC/"), varnameP  = "pr_sample",
                       ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r14_hist_ReGrid/"),
                                    paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r14_hist_ReGrid/")))
            ,
            EOC=list(name="EOC",col="red", stam = "his.CXFPS025.HCLIM38h1_CXFPS_fRACMOfECEARTH_r04_rcp85.AYEAR.1hr.G5x5",
                     maskdir = paste0(hdir,"/analysis/CXFPS/PREPROCESS/NC/"), varnameP  = "pr_mean",
                     ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r04_rcp85_ReGrid/"),
                                  paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r04_rcp85_ReGrid/")))
            ,
            EOC2=list(name="EOC",col="red", stam = "his.CXFPS025.HCLIM38h1_CXFPS_fRACMOfECEARTH_r04_rcp85.AYEAR.1hr.G5x5",
                      maskdir = paste0(hdir,"/analysis/CXFPS/PREPROCESS/NC/"), varnameP  = "pr_sample",
                      ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r04_rcp85_ReGrid/"),
                                   paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS_fRACMOfECEARTH_r04_rcp85_ReGrid/")))

          )
          varnameT  <- "td2m_mean"
          varnameTas  <- "tas_mean"
        }

        if(PROJDOM=="CPM-ERAINT-HCLIM-NWE3"){
          EXPLIST = list(
            NWE3CTL=list(name="CTL",col="black", stam = "his.EUCPNW025.HCLIM38h1_EUCP_CTL_fRACMOfERAIN.AYEAR.1hr.G5x5",
                         maskdir = paste0(hdir,"/analysis/EUCP_NW/PREPROCESS/NC/"), varnameP  = "pr_mean",
                         ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/EUCP_NW_REGRID/HCLIM38h1_EUCP_CTL_fRACMOfERAIN_ReGrid/"),
                                      paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/EUCP_NW_REGRID/HCLIM38h1_EUCP_CTL_fRACMOfERAIN_ReGrid/")),
                         yearsel = c(2008,2018))
            ,
            NWE3CTL2=list(name="CTL",col="black", stam = "his.EUCPNW025.HCLIM38h1_EUCP_CTL_fRACMOfERAIN.AYEAR.1hr.G5x5",
                          maskdir = paste0(hdir,"/analysis/EUCP_NW/PREPROCESS/NC/"), varnameP  = "pr_sample",
                          ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/EUCP_NW_REGRID/HCLIM38h1_EUCP_CTL_fRACMOfERAIN_ReGrid/"),
                                       paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/EUCP_NW_REGRID/HCLIM38h1_EUCP_CTL_fRACMOfERAIN_ReGrid/")),
                          yearsel = c(2008,2018))
            ,
            NWE3PGW=list(name="PGW",col="red", stam = "his.EUCPNW025.HCLIM38h1_EUCP_DTglob1H_ECEARTH_fRACMOfERAIN.AYEAR.1hr.G5x5",
                         maskdir = paste0(hdir,"/analysis/EUCP_NW/PREPROCESS/NC/"), varnameP  = "pr_mean",
                         ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/EUCP_NW_REGRID/HCLIM38h1_EUCP_DTglob1H_ECEARTH_fRACMOfERAIN_ReGrid/"),
                                      paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX//EUCP_NW_REGRID/HCLIM38h1_EUCP_DTglob1H_ECEARTH_fRACMOfERAIN_ReGrid/")),
                         yearsel = c(2008,2018))
            ,
            NWE3PGW2=list(name="PGW",col="red", stam = "his.EUCPNW025.HCLIM38h1_EUCP_DTglob1H_ECEARTH_fRACMOfERAIN.AYEAR.1hr.G5x5",
                           maskdir = paste0(hdir,"/analysis/EUCP_NW/PREPROCESS/NC/"), varnameP  = "pr_sample",
                          ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/EUCP_NW_REGRID/HCLIM38h1_EUCP_DTglob1H_ECEARTH_fRACMOfERAIN_ReGrid/"),
                                       paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX//EUCP_NW_REGRID/HCLIM38h1_EUCP_DTglob1H_ECEARTH_fRACMOfERAIN_ReGrid/")),
                          yearsel = c(2008,2018))

          )
          varnameT    <- "td2m_mean"
          varnameTas  <- "tas_mean"
        }

        if(PROJDOM=="CPM-ERAINT-HCLIM-ALP3"){
          EXPLIST = list(
            NWE3CTL=list(name="CTL",col="black", stam = "his.CXFPS025.HCLIM38h1_CXFPS1999_fRACMO.AYEAR.1hr.G5x5",
                         maskdir = paste0(hdir,"/analysis/CXFPS/PREPROCESS/NC/"), varnameP  = "pr_mean",
                         ncindir =  c(paste0(hdir,"/LINKS/PC132057_NOB2/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS1999_fRACMO_ReGrid/"),
                                      paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS1999_fRACMO_ReGrid/"))),
            NWE3CTL2=list(name="CTL",col="black", stam = "his.CXFPS025.HCLIM38h1_CXFPS1999_fRACMO.AYEAR.1hr.G5x5",
                          maskdir = paste0(hdir,"/analysis/CXFPS/PREPROCESS/NC/"), varnameP  = "pr_sample",
                          ncindir =  c(paste0(hdir,"/LINKS/PC132057_NOB2/CORDEX/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS1999_fRACMO_ReGrid/"),
                                       paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX_FPS_REGRID/HCLIM38h1_CXFPS1999_fRACMO_ReGrid/")))

          )
          varnameT  <- "td2m_mean"
          varnameTas  <- "tas_mean"
        }

        if(PROJDOM=="CPM-ERAINT-ETH-ALP3"){
          EXPLIST = list(
            CTL=list(name="CTL",col="black", stam = "ALP-3_ECMWF-ERAINT_evaluation_r1i1p1_COSMO-pompa_5.0_2019.1_1hr_2000_2009_remapGL.AYEAR.2000-2009.G5x5",
                     maskdir = paste0(hdir,"/analysis/CXFPS_ETH/PREPROCESS/NC/"), varnameP  = "pr_mean",
                     ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/ETH_CXFPS_ERA_Regrid/"),
                                  paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/ETH_CXFPS_ERA_Regrid/"))),
            CTL2=list(name="CTL",col="black", stam = "ALP-3_ECMWF-ERAINT_evaluation_r1i1p1_COSMO-pompa_5.0_2019.1_1hr_2000_2009_remapGL.AYEAR.2000-2009.G5x5",
                      maskdir = paste0(hdir,"/analysis/CXFPS_ETH/PREPROCESS/NC/"), varnameP  = "pr_sample",
                      ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/ETH_CXFPS_ERA_Regrid/"),
                                   paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/ETH_CXFPS_ERA_Regrid/")))

          )
          varnameT  <- "td2m_mean"
          varnameTas  <- "tas_mean"
          pmul = 1.
        }

        if(PROJDOM=="CPM-ERAINT-METO-ALP3"){
          EXPLIST = list(
            CTL=list(name="CTL",col="black", stam = "CER-2p2_ECMWF-ERAINT_evaluation_r1i1p1_UKMO-UM10p1_v01.AYEAR.2000-2011.G5x5",
                     maskdir = paste0(hdir,"/analysis/CXFPS_METO/PREPROCESS/NC/"), varnameP  = "pr_mean",
                     ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/METO_CXFPS_ERA_ReGrid/"),
                                  paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/METO_CXFPS_ERA_ReGrid/"))),
            CTL2=list(name="CTL",col="black", stam = "CER-2p2_ECMWF-ERAINT_evaluation_r1i1p1_UKMO-UM10p1_v01.AYEAR.2000-2011.G5x5",
                      maskdir = paste0(hdir,"/analysis/CXFPS_METO/PREPROCESS/NC/"), varnameP  = "pr_sample",
                      ncindir =  c(paste0(hdir,"/LINKS/PC190611_NOB/CORDEX/CORDEX_FPS_REGRID/METO_CXFPS_ERA_ReGrid/"),
                                   paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/METO_CXFPS_ERA_ReGrid/")))

          )
          varnameT  <- "td2m_mean"
          varnameTas  <- "tas_mean"
        }


        if (PROJDOM=="CPM-ERAINT-CNRM-NWE3"){

          EXPLIST = list(
            CTL=list(name="CTL",col="black", stam = "his.NWE-3_ECMWF-ERAINT_evaluation_r1i1p1_CNRM-AROME41t1_fpsconv-x2yn2-v1_1hr.AYEAR.1hr.G5x5",
                     maskdir = paste0(hdir,"/analysis/CXFPS_CNRM/PREPROCESS/NC/"), varnameP  = "pr_mean",
                     ncindir =  c(paste0(hdir,"/LINKS/ddisk2/CORDEX/CORDEX_FPS_REGRID/CNRM_ERA_ReGrid/"),
                                  paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/CNRM_ERA_ReGrid//"))),
            CTL2=list(name="CTL",col="black", stam = "his.NWE-3_ECMWF-ERAINT_evaluation_r1i1p1_CNRM-AROME41t1_fpsconv-x2yn2-v1_1hr.AYEAR.1hr.G5x5",
                      maskdir = paste0(hdir,"/analysis/CXFPS_CNRM/PREPROCESS/NC/"), varnameP  = "pr_sample",
                      ncindir =  c(paste0(hdir,"/LINKS/ddisk2/CORDEX/CORDEX_FPS_REGRID/CNRM_ERA_ReGrid/"),
                                   paste0(hdir,"/LINKS/DATA_EXT_ROODWIT/CORDEX/CORDEX_FPS_REGRID/CNRM_ERA_ReGrid//")))

          )
          varnameT  <- "td2m_mean"
          varnameTas  <- "tas_mean"
        }






        # ------------------------------------------------------------------------------------------------------------------------------
        # start analysis
        # ------------------------------------------------------------------------------------------------------------------------------

        for (iperiod in 1:length(EXPLIST)) {


          varnameP  = paste0(EXPLIST[[iperiod]]$varnameP)

          # read in model/obs data

          nquant <- length(pquant)
          
          ncindir = "./"
          ltest <- dir.exists(EXPLIST[[iperiod]]$ncindir)
          if (any(ltest)) {
            ncindir <- EXPLIST[[iperiod]]$ncindir[min(which(ltest))]
          } else {
            print("no valid directory for netcdf files found; WARNING ; can be neglected here !!")
          }

          maskfile <- paste0(ncindir,EXPLIST[[iperiod]]$maskdir,"/mask_",maskid,".nc") # first relative path
          if ( ! file.exists(maskfile) )  {maskfile <- paste0(EXPLIST[[iperiod]]$maskdir,"/mask_",maskid,".nc")} # try absolut path

#          if ( ! file.exists(maskfile) & ! maskfile == "no" ) {
           if ( ! file.exists(maskfile) & ! maskfile == "no" & lRdatreprocess  ) {
              
              print(paste("maskfile does not exist",maskfile))

          } else {

            ltest <- dir.exists(EXPLIST[[iperiod]]$ncindir)
            if (any(ltest)) {
              ncindir <- EXPLIST[[iperiod]]$ncindir[min(which(ltest))]
            } else {
              print("no valid directory for netcdf files found; WARNING, but can be neglected here !!")
            }

              tas_infile <- paste0(ncindir,paste0(varnameTas,"_",EXPLIST[[iperiod]]$stam,"_3hr_",EXPLIST[[iperiod]]$name,".nc"))
              td_infile <-  paste0(ncindir,paste0(varnameT,"_",EXPLIST[[iperiod]]$stam,"_3hr_",EXPLIST[[iperiod]]$name,".nc"))
              pr_infile <-  paste0(ncindir,paste0(varnameP,"_",EXPLIST[[iperiod]]$stam,"_1hr_",EXPLIST[[iperiod]]$name,".nc"))

              # tries to read Rdata file if exists
              RdataFile <- paste0(ncindir,"/PrTdTas_",maskid,"_",varnameP,"_",EXPLIST[[iperiod]]$name,".Rdata")
              
              # this part is a shortcut to the saved and extracted Rdata files
              rddir = paste0(rodir,PROJDOM,"/")
              RdataFile <- paste0(rddir,"/PrTdTas_",maskid,"_",varnameP,"_",EXPLIST[[iperiod]]$name,".Rdata")
              
 
              if ( ! file.exists(RdataFile) | lRdatreprocess ){

                  indataP <- gl.ReadFromNetcdf2D_Tsteps(ncinfile = pr_infile,varnameP,maskfile = maskfile,Nblock = 8)
                  indataT <- gl.ReadFromNetcdf2D_Tsteps(ncinfile = td_infile,varnameT,maskfile = maskfile,Nblock = 8)
                  indataTas <- gl.ReadFromNetcdf2D_Tsteps(ncinfile = tas_infile,varnameTas,maskfile = maskfile,Nblock = 8)
                  save(indataP,indataT,indataTas,file=RdataFile)

              } else {

                load(RdataFile)
              }

            # read in data
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
              pr_season = pr_season * pmul
              pr_season[pr_season > -0.05 & pr_season <0.] <- 0. # some model output contain negative because of acc, rain stored

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
              nbsize = 10 # sample size of bootstrap resample
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

          plotfile <-  paste0(dirplots,"/Cluster0_AllStatBootstrap_",dataid,".eps")
          plotfile <-  paste0(dirplots,"/Cluster0_AllStatBootstrap_",dataid)
          gl.PlotScalingBoot(ScalingList,plotfile=plotfile)

        } # iperiod



      } # end loop of project


    } # seasonid
} # maskid
