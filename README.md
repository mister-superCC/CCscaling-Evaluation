# CCscaling-Evaluation
Scripts needed for Lenderink et al 2024

ScalingLibSmall.R contains a subset of my routines to create scaling object, statistics and reading and plotting of data

scaling for observations
ScalingCirculationIndexStandAlone_obs.R 

The data for NL are freely available
The data for SFR are NOT freely available. The processed object files (.Rdata) which are produced in directories data_hourly_bootstrap data_hourly_bootstrap_abs can be obtained here: ...

scaling for models
ScalingCirculationIndexV3_bootstrap_models.R
Runs from data here: 
This data has been already processed to contain only the grid points associated with the stations studied. Scripts can be run from the regridded model data.

Conversion new output format object to old object
Run with pr_sample and pr_mean
conversion_boot2noboot.R

Make most of the plots
MergeResultsBootstrap.R

Some other plots made by:
CompareFromBootstrap_SFRcomparedNL.R
ComparePOE_NL2SFR.R

MergePlotsFromScaling_bootstrap_MJJAS.sc
contains the merge of the different subplots to the plots in the paper. Uses a simple Latex based program to merge the subplots (not provided as it contains some systems specific binaries) but you can see here which are the subplots
