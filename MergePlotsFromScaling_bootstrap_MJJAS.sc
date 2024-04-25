#!/bin/csh -f


set opt = 1  

set lexit = F
set lepsrm = T
set season = SHY
set season = MJJAS
# set season = AMJJ # ASON # AMJJ 
set cpmprname = pr_sample
#set cpmprname = pr_mean
set indir  = plots_merged_bootstrap/
set indir2 = plots_scaling/
set indir3 = plots_CDF_NLvsSFR/

set plotid = $1

# SFRstations.land.all SFRstationsMedSea.land.all SFRstationsNoMedSea.land.all NLstations.land.all

if ( $opt == 1) then
  set odir = $indir/catted_threeareas/
 set pctl1 = P99
 set pctl2 = P99.9
 set areas = ( NLstations.land.all  SFRstationsNoMedSea.land.all SFRstationsMedSea.land.all  )
endif



if ( ! -r $odir) mkdir $odir
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/
if ( ! -r $odir) mkdir $odir



goto $plotid

# paper plot 2, CFD extremes of all areas

plot2:
set plotid = plot2
# proobability of exceedence plot
set area1 = $areas[1]
set area2 = $areas[2]
set area3 = $areas[3]
set ptype = CDF

#set area = simpleNL.all.oroglt400

set pfiles = ( \
               $indir/CPM_${cpmprname}_${season}-${area1}_${ptype}.eps \
               $indir/CPM_${cpmprname}_${season}-${area2}_${ptype}.eps \
               $indir/CPM_${cpmprname}_${season}-${area3}_${ptype}.eps \
               $indir/RCM_pr_${season}-${area1}_${ptype}.eps \
               $indir/RCM_pr_${season}-${area2}_${ptype}.eps \
               $indir/RCM_pr_${season}-${area3}_${ptype}.eps \
               )
                
               ls -l $pfiles
                
                
set ofile = $odir/PoE6P_${season}_${ptype}.$plotid.eps
multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
epstopdf $ofile
if ($lepsrm == "T") rm $ofile

exit

 
 
#paper plot3

# 50,80th percentile of distribution of TD and DPD as function of intensity

plot3:
set plotid = plot3
foreach area ( $areas ) 

		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_CDFDPdep_P50_TD-${area}_ptype2.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFDPdep_P50-${area}_ptype1.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFDPdep_P80-${area}_ptype1.eps \
					   $indir/RCM_pr_${season}_CDFDPdep_P50_TD-${area}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFDPdep_P50-${area}_ptype1.eps \
					   $indir/RCM_pr_${season}_CDFDPdep_P80-${area}_ptype1.eps \
									 )
		 
		set ofile = $odir/CDFDPdep_${season}_${area}.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile

end



exit


plot4:
# full CDF per statistics, 

# RCM_pr_SON_CDFsingle_P99_anom-SFRstationsMedSea.land.all_ptype2.eps

set plotid = plot4
set stats = ( "P99_anom"  "all_anom" "wet_anom" )

foreach stat ( $stats ) 

		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stat}-${areas[1]}_ptype2.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stat}-${areas[2]}_ptype2.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stat}-${areas[3]}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stat}-${areas[1]}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stat}-${areas[2]}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stat}-${areas[3]}_ptype2.eps \
									 )
		 
		ls -l $pfiles
	
		set ofile = $odir/CDFsingle_${season}_$stat.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile

end

exit


plot4_sup:
# full CDF per statistics, 

# RCM_pr_SON_CDFsingle_P99_anom-SFRstationsMedSea.land.all_ptype2.eps

set plotid = plot4
set stats = ( "P99_anom_TD"  "all_anom_TD" "wet_anom_TD" )

foreach stat ( $stats ) 

		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stat}-${areas[1]}_ptype4.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stat}-${areas[2]}_ptype4.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stat}-${areas[3]}_ptype4.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stat}-${areas[1]}_ptype4.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stat}-${areas[2]}_ptype4.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stat}-${areas[3]}_ptype4.eps \
									 )
		 
		ls -l $pfiles
	
		set ofile = $odir/CDFsingle_${season}_$stat.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile

end

exit



# plot5, scaling diagram based on observations



plot5:

# plots_scaling/ScalingTD_precip-OBS-MEDSEA-STATIONS-MEDSEAstations-SFRstationsMedSea.land.all-SHY-1991-2020.eps

set plotid = plot5
set ptype = ScalingTD

set pfiles = ( \
   $indir2/${ptype}_precip-OBS-NL-STATIONS-NLstations-NLstations.land.all-${season}-1991-2020.eps \
   $indir2/${ptype}_precip-OBS-NoMEDSEA-STATIONS-NoMEDSEAstations-SFRstationsNoMedSea.land.all-${season}-1991-2020.eps \
   $indir2/${ptype}_precip-OBS-MEDSEA-STATIONS-MEDSEAstations-SFRstationsMedSea.land.all-${season}-1991-2020.eps \
)

ls -l $pfiles
set ofile = $odir/${ptype}_${season}_OBS.$plotid.eps
multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
epstopdf $ofile
echo $ofile
if ($lepsrm == "T") rm $ofile


exit




# TD scaling conidtional on DPD

plot6:

set plotid = plot6

set ptype1 = ptype6_whf_lowRH
set ptype2 = ptype6_whf_highRH

#RCM_pr_SHY_P99-SFRstationsNoMedSea.land.all_ptype1_RHlow.eps

foreach area ( $areas )
		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_P99-${area}_ptype1_RHall.eps \
					   $indir/CPM_${cpmprname}_${season}_P99-${area}_ptype1_RHhigh.eps \
					   $indir/CPM_${cpmprname}_${season}_P99-${area}_ptype1_RHlow.eps \
					   $indir/RCM_pr_${season}_P99-${area}_ptype1_RHall.eps \
					   $indir/RCM_pr_${season}_P99-${area}_ptype1_RHhigh.eps \
					   $indir/RCM_pr_${season}_P99-${area}_ptype1_RHlow.eps \
									 )
		 
		set ofile = $odir/ScalingRHdep_${season}_${area}.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
        echo $ofile
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile
end


set ptype1 = ptype6_whf_lowRH
set ptype2 = ptype6_whf_highRH

#RCM_pr_SHY_P99-SFRstationsNoMedSea.land.all_ptype1_RHlow.eps

foreach area ( $areas )
		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_P99.9-${area}_ptype1_RHall.eps \
					   $indir/CPM_${cpmprname}_${season}_P99.9-${area}_ptype1_RHhigh.eps \
					   $indir/CPM_${cpmprname}_${season}_P99.9-${area}_ptype1_RHlow.eps \
					   $indir/RCM_pr_${season}_P99.9-${area}_ptype1_RHall.eps \
					   $indir/RCM_pr_${season}_P99.9-${area}_ptype1_RHhigh.eps \
					   $indir/RCM_pr_${season}_P99.9-${area}_ptype1_RHlow.eps \
									 )
		 
		set ofile = $odir/ScalingRHdep_${season}_P999_${area}.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
        echo $ofile
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile
end


exit




exit


plot7:
# scatter plots Pat15c versus Scaling

set plotid = plot7
set pctl = P99
set RHs = ( allRH lowRH highRH  )
# set season = SHY

# need to pull in season

set allfiles = ()
foreach RH ( $RHs )
    set pfiles = ( $indir2/AlphaVSPat15_NL_${RH}_${pctl}_$season.eps \
                   $indir2/AlphaVSPat15_SFR-cent_${RH}_${pctl}_$season.eps \
		           $indir2/AlphaVSPat15_SFR-med_${RH}_${pctl}_$season.eps )
	set allfiles = ( $allfiles $pfiles )

#	epstopdf AlphaVSPat15_AllAreas_${RH}_${pctl}_$season.eps
end

set ofile = $odir/Pat15vsScaling_${season}_${pctl}.$plotid.eps
multeps2one_ls.sc -f $allfiles -xy 3 3 -w 5 -eps -o $ofile -clean  
echo $ofile

epstopdf $ofile
if ($lepsrm == "T") rm $ofile

# versie 2
set plotid = plot7
set pctl = P99
set RHs = ( allRH lowRH highRH  )
# set season = SHY

# need to pull in season

set allfiles = ()
foreach areaid ( NL SFR-cent )
    set pfiles = ( $indir2/AlphaVSPat15_${areaid}_allRH_${pctl}_$season.eps \
			       $indir2/AlphaVSPat15_${areaid}_highRH_${pctl}_$season.eps  \
			       $indir2/AlphaVSPat15_${areaid}_lowRH_${pctl}_$season.eps )
	set allfiles = ( $allfiles $pfiles )
end

set ofile = $odir/Pat15vsScaling_${season}_${pctl}.${plotid}_v2.eps
multeps2one_ls.sc -f $allfiles -xy 3 2 -w 5 -eps -o $ofile -clean  
echo $ofile

epstopdf $ofile
if ($lepsrm == "T") rm $ofile


exit



# plot intensity low RH versus high RH as function of dew point


plot8:
set plotid = plot8

foreach area ( $areas ) 
		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_P99-${area}_ptype2_RHlow_vs_high.eps \
					   $indir/RCM_pr_${season}_P99-${area}_ptype2_RHlow_vs_high.eps \
									  )
		 
		set ofile = $odir/ScalingRHlowVShigh_${season}_${area}.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 2 1 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile
end

exit



# paper plot wet hour fraction



plot9:


set plotid = plot9
set ptype3 = ptype6_whf_lowRH
set ptype2 = ptype6_whf_highRH
set ptype1 = ptype6_whf_all


foreach area ( $areas ) 

		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}-${area}_${ptype1}.eps \
					   $indir/CPM_${cpmprname}_${season}-${area}_${ptype2}.eps \
 					   $indir/CPM_${cpmprname}_${season}-${area}_${ptype3}.eps \
                       $indir/RCM_pr_${season}-${area}_${ptype1}.eps \
					   $indir/RCM_pr_${season}-${area}_${ptype2}.eps \
			  	       $indir/RCM_pr_${season}-${area}_${ptype3}.eps \
									  )
						
					   ls -l $pfiles
						
						
		set ofile = $odir/WHF_${season}_${area}.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile
		
end

exit


plot10:

set plotid = plot10

set pfiles = ( \
	$indir3/CPM_pr_sample_NoMEDSEAvsNL_PoE_$season.eps \
 	$indir3/RCM_pr_NoMEDSEAvsNL_PoE_$season.eps )
#       $indir3/CPM_pr_sample_MEDSEAvsNL_PoE_$season.eps \
#   	$indir3/RCM_pr_MEDSEAvsNL_PoE_$season.eps )
	
	
set ofile = $odir/POE_comparison_${season}_SFRsNL.$plotid.eps
multeps2one_ls.sc -f $pfiles -xy 2 1 -w 8 -eps -o $ofile -clean  
epstopdf $ofile
if ($lepsrm == "T") rm $ofile


exit		
		
	
# ==================================================================================================================================================================
# seperate plots for the supplement; note some of the plots for the supplement are already produced by the normal plotting above, e.g. Figure 1 for NL and SFR-med




plot2_sup:
# full CDF per statistics, 

# RCM_pr_SON_CDFsingle_P99_anom-SFRstationsMedSea.land.all_ptype2.eps

set stats = ( "all_anom_TD" "wet_anom_TD"  "P99_anom_TD"  )

foreach area ( $areas ) 

		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stats[1]}-${area}_ptype4.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stats[2]}-${area}_ptype4.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stats[3]}-${area}_ptype4.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stats[1]}-${area}_ptype4.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stats[2]}-${area}_ptype4.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stats[3]}-${area}_ptype4.eps \
									 )
		 
		ls -l $pfiles
	
		set ofile = $odir/CDFsingle_${season}_$area.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile

end

exit


plot3_sup:
# full CDF per statistics, 

# RCM_pr_SON_CDFsingle_P99_anom-SFRstationsMedSea.land.all_ptype2.eps

set stats = ( "all_anom" "wet_anom"  "P99_anom"  )

foreach area ( $areas ) 

		set pfiles = ( \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stats[1]}-${area}_ptype2.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stats[2]}-${area}_ptype2.eps \
					   $indir/CPM_${cpmprname}_${season}_CDFsingle_${stats[3]}-${area}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stats[1]}-${area}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stats[2]}-${area}_ptype2.eps \
					   $indir/RCM_pr_${season}_CDFsingle_${stats[3]}-${area}_ptype2.eps \
									 )
		 
		ls -l $pfiles
	
		set ofile = $odir/CDFsingle_${season}_$area.$plotid.eps
		multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
		epstopdf $ofile
		if ($lepsrm == "T") rm $ofile

end

exit





plot5_abs:

# plots_scaling/ScalingTD_precip-OBS-MEDSEA-STATIONS-MEDSEAstations-SFRstationsMedSea.land.all-SHY-1991-2020.eps


set indir  = plots_merged_bootstrap_abs/
set indir2 = plots_scaling_abs/
set odir = $indir/catted_threeareas/


if ( ! -r $odir) mkdir $odir
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/
if ( ! -r $odir) mkdir $odir


set ptype = ScalingTD

set pfiles = ( \
   $indir2/${ptype}_precip-OBS-NL-STATIONS-NLstations-NLstations.land.all-${season}-1991-2020.eps \
   $indir2/${ptype}_precip-OBS-NoMEDSEA-STATIONS-NoMEDSEAstations-SFRstationsNoMedSea.land.all-${season}-1991-2020.eps \
   $indir2/${ptype}_precip-OBS-MEDSEA-STATIONS-MEDSEAstations-SFRstationsMedSea.land.all-${season}-1991-2020.eps \
)

ls -l $pfiles
set ofile = $odir/${ptype}_${season}_OBS.$plotid.eps
multeps2one_ls.sc -f $pfiles -xy 3 2 -w 8 -eps -o $ofile -clean  
epstopdf $ofile
echo $ofile
if ($lepsrm == "T") rm $ofile


exit


plot7_abs:


# versie 2
set pctl = P99.9
set indir  = plots_merged_bootstrap_abs/
set indir2 = plots_scaling_abs/

set RHs = ( allRH lowRH highRH  )
# set season = SHY

# need to pull in season
set odir = $indir/catted_threeareas/
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/

set allfiles = ()
foreach areaid ( NL SFR-cent )
    set pfiles = ( $indir2/AlphaVSPat15_${areaid}_allRH_${pctl}_$season.eps \
			       $indir2/AlphaVSPat15_${areaid}_highRH_${pctl}_$season.eps  \
			       $indir2/AlphaVSPat15_${areaid}_lowRH_${pctl}_$season.eps )
	set allfiles = ( $allfiles $pfiles )
end

set ofile = $odir/Pat15vsScaling_${season}_${pctl}.${plotid}_v2.eps
multeps2one_ls.sc -f $allfiles -xy 3 2 -w 5 -eps -o $ofile -clean  
echo $ofile

epstopdf $ofile
if ($lepsrm == "T") rm $ofile


exit


plot5_sup_RCMs:

set indir  = plots_merged_bootstrap/
set indir2 = plots_scaling/
set odir = $indir/catted_threeareas/

if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/


set ifiles = ( \
$indir2/ScalingTD_precip-OBS-NL-STATIONS-NLstations-NLstations.land.all-${season}-1991-2020.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-CLMcom-ETH-COSMO-NLstations.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-CNRM-ALADIN63-NLstations.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-DMI-HIRHAM5-NLstations.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-GERICS-REMO2015-NLstations.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-KNMI-RACMO22E-NLstations.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-MOHC-HadREM3-NLstations.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-SMHI-RCA4-NLstations.land.all-${season}-1991-2010.eps \
)

set ofile = $odir/ScalingRCMs_NL_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 3 3 -w 5 -eps -o $ofile -clean  
echo $ofile

epstopdf $ofile
if ($lepsrm == "T") rm $ofile


set ifiles = ( \
$indir2/ScalingTD_precip-OBS-NoMEDSEA-STATIONS-NoMEDSEAstations-SFRstationsNoMedSea.land.all-${season}-1991-2020.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-CLMcom-ETH-COSMO-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-CNRM-ALADIN63-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-DMI-HIRHAM5-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-GERICS-REMO2015-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-KNMI-RACMO22E-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-MOHC-HadREM3-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
$indir2/ScalingTD_pr-RCM-ERAINT-SMHI-RCA4-SFRstationsNoMedSea.land.all-${season}-1991-2010.eps \
)

set ofile = $odir/ScalingRCMs_SFR-cent_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 3 3 -w 5 -eps -o $ofile -clean  
echo $ofile

epstopdf $ofile
if ($lepsrm == "T") rm $ofile

exit

plot5_sup_CPMs:


set indir  = plots_merged_bootstrap/
set indir2 = plots_scaling/
set odir = $indir/catted_threeareas/

if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/


set ifiles = ( \
$indir2/ScalingTD_precip-OBS-NL-STATIONS-NLstations-NLstations.land.all-${season}-1991-2020.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-CNRM-NWE3-NLstations.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-ETH-ALP3-NLstations.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-HCLIM-ALP3-NLstations.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-HCLIM-NWE3-NLstations.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-METO-ALP3-NLstations.land.all-${season}-CTL.eps \
)

set ofile = $odir/ScalingCPMs_NL_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 3 2 -w 5 -eps -o $ofile -clean  
echo $ofile
epstopdf $ofile
if ($lepsrm == "T") rm $ofile



set ifiles = ( \
$indir2/ScalingTD_precip-OBS-NoMEDSEA-STATIONS-NoMEDSEAstations-SFRstationsNoMedSea.land.all-${season}-1991-2020.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-CNRM-NWE3-SFRstationsNoMedSea.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-ETH-ALP3-SFRstationsNoMedSea.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-HCLIM-ALP3-SFRstationsNoMedSea.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-HCLIM-NWE3-SFRstationsNoMedSea.land.all-${season}-CTL.eps \
$indir2/ScalingTD_${cpmprname}-CPM-ERAINT-METO-ALP3-SFRstationsNoMedSea.land.all-${season}-CTL.eps \
)

set ofile = $odir/ScalingCPMs_SFR-cent_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 3 2 -w 5 -eps -o $ofile -clean  
echo $ofile
epstopdf $ofile
if ($lepsrm == "T") rm $ofile




plot13_sup_humdist:

set indir  = plots_merged_bootstrap/
set odir = $indir/catted_threeareas/
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/

set dataid = OBS


set ifiles = ( \
$indir/${dataid}_CDFDPsingle_allstat_TD_precip_${season}-NLstations.land.all_ptype5.eps    \
$indir/${dataid}_CDFDPsingle_allstat_TD_precip_${season}-SFRstationsNoMedSea.land.all_ptype5.eps \
$indir/${dataid}_CDFDPsingle_allstat_precip_${season}-NLstations.land.all_ptype3.eps \
$indir/${dataid}_CDFDPsingle_allstat_precip_${season}-SFRstationsNoMedSea.land.all_ptype3.eps \
)

set ofile = $odir/Distribution_${dataid}_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 2 2 -w 5 -eps -o $ofile -clean  
echo $ofile
epstopdf $ofile
if ($lepsrm == "T") rm $ofile

#exit


foreach dataid ( ALADIN CLMCOM HADRM HIRHAM RACMO RCAO REMO ) 

    set ifiles = ( \
    $indir/${dataid}_CDFDPsingle_allstat_TD_pr_${season}-NLstations.land.all_ptype5.eps    \
    $indir/${dataid}_CDFDPsingle_allstat_TD_pr_${season}-SFRstationsNoMedSea.land.all_ptype5.eps \
    $indir/${dataid}_CDFDPsingle_allstat_pr_${season}-NLstations.land.all_ptype3.eps \
    $indir/${dataid}_CDFDPsingle_allstat_pr_${season}-SFRstationsNoMedSea.land.all_ptype3.eps \
    )

    set ofile = $odir/Distribution_${dataid}_${season}_${plotid}.eps
    multeps2one_ls.sc -f $ifiles -xy 2 2 -w 5 -eps -o $ofile -clean  
    echo $ofile
    epstopdf $ofile
    if ($lepsrm == "T") rm $ofile

end



foreach dataid ( ALADIN CLMCOM HADRM HIRHAM RACMO RCAO REMO ) 

    set ifiles = ( \
    $indir/${dataid}_CDFDPsingle_allstat_TD_pr_${season}-NLstations.land.all_ptype5.eps    \
    $indir/${dataid}_CDFDPsingle_allstat_TD_pr_${season}-SFRstationsNoMedSea.land.all_ptype5.eps \
    $indir/${dataid}_CDFDPsingle_allstat_pr_${season}-NLstations.land.all_ptype3.eps \
    $indir/${dataid}_CDFDPsingle_allstat_pr_${season}-SFRstationsNoMedSea.land.all_ptype3.eps \
    )

    set ofile = $odir/Distribution_${dataid}_${season}_${plotid}.eps
    multeps2one_ls.sc -f $ifiles -xy 2 2 -w 5 -eps -o $ofile -clean  
    echo $ofile
    epstopdf $ofile
    if ($lepsrm == "T") rm $ofile

end


exit


plot15_abs_scalingRHdep:

set indir2 =  plots_scaling_abs/

   set ifiles = ( \
    $indir2/CPM_${cpmprname}_${season}_P99.9-SFRstationsNoMedSea.land.all_ptype1_RHhigh.eps    \
    $indir2/CPM_${cpmprname}_${season}_P99.9-SFRstationsNoMedSea.land.all_ptype1_RHlow.eps \
    $indir2/RCM_pr_${season}_P99.9-SFRstationsNoMedSea.land.all_ptype1_RHhigh.eps    \
    $indir2/RCM_pr_${season}_P99.9-SFRstationsNoMedSea.land.all_ptype1_RHlow.eps \
   )



set indir  = plots_merged_bootstrap_abs/
set odir = $indir/catted_threeareas/
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/

set ofile = $odir/Scaling_SFR-cent_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 2 2 -w 5 -eps -o $ofile -clean  
echo $ofile
epstopdf $ofile
if ($lepsrm == "T") rm $ofile



  set ifiles = ( \
    $indir2/CPM_${cpmprname}_${season}_P99.9-NLstations.land.all_ptype1_RHhigh.eps    \
    $indir2/CPM_${cpmprname}_${season}_P99.9-NLstations.land.all_ptype1_RHlow.eps \
    $indir2/RCM_pr_${season}_P99.9-NLstations.land.all_ptype1_RHhigh.eps    \
    $indir2/RCM_pr_${season}_P99.9-NLstations.land.all_ptype1_RHlow.eps \
   )



set indir  = plots_merged_bootstrap_abs/
set odir = $indir/catted_threeareas/
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/

set ofile = $odir/Scaling_NL_${season}_${plotid}.eps
multeps2one_ls.sc -f $ifiles -xy 2 2 -w 5 -eps -o $ofile -clean  
echo $ofile
epstopdf $ofile
if ($lepsrm == "T") rm $ofile




set indir2 =  plots_scaling_abs/

   set ifiles = ( \
    $indir2/CPM_${cpmprname}_${season}_P99.99-SFRstationsNoMedSea.land.all_ptype1_RHhigh.eps    \
    $indir2/CPM_${cpmprname}_${season}_P99.99-SFRstationsNoMedSea.land.all_ptype1_RHlow.eps \
    $indir2/RCM_pr_${season}_P99.99-SFRstationsNoMedSea.land.all_ptype1_RHhigh.eps    \
    $indir2/RCM_pr_${season}_P99.99-SFRstationsNoMedSea.land.all_ptype1_RHlow.eps \
   )



set indir  = plots_merged_bootstrap_abs/
set odir = $indir/catted_threeareas/
if ( $cpmprname == "pr_mean" )   set odir = $odir/pr_mean/
if ( $cpmprname == "pr_sample" ) set odir = $odir/pr_sample/

set ofile = $odir/Scaling_SFR-cent_${season}_${plotid}_P99.99.eps
multeps2one_ls.sc -f $ifiles -xy 2 2 -w 5 -eps -o $ofile -clean  
echo $ofile
epstopdf $ofile
if ($lepsrm == "T") rm $ofile







