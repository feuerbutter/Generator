#------------------------------------------------------------------------------------------
# Submit a GENIE/SK event generation job using a histogram-based neutrino flux description
#
# Syntax:
#   shell% perl submit-evg_sk_fhst.pl <options>
#
# Options:
#    --version         : GENIE version
#    --run             : 1,2,3,...
#    --neutrino        : numu, numubar, nue, nuesig
#   [--flux-version]   : JNUBEAM flux version, <07a, 10>, default: 10
#   [--flux-hist-file] : JNUBEAM flux histogram file, default: sk_flux_histograms.root
#   [--arch]           : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]     : default: <version>
#   [--cycle]          : default: 01
#   [--use-valgrind]   : default: off
#   [--batch-system]   : <PBS, LSF>, default: PBS
#   [--queue]          : default: prod
#   [--softw-topdir]   : default: /opt/ppd/t2k/GENIE
#
# Example:
#   shell& perl submit-evg_sk_fhst.pl --run 180 --neutrino numubar --version v2.5.1
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#------------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version  = $ARGV[$iarg+1]; }
  if($_ eq '--run'     )       { $run            = $ARGV[$iarg+1]; }
  if($_ eq '--neutrino')       { $neutrino       = $ARGV[$iarg+1]; }
  if($_ eq '--flux-version')   { $flux_version   = $ARGV[$iarg+1]; }
  if($_ eq '--flux-hist-file') { $flux_hist_file = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch           = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production     = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle          = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind   = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system   = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue          = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir   = $ARGV[$iarg+1]; }  
  $iarg++;
}
die("** Aborting [Undefined run number. Use the --run option]")
unless defined $run;
die("** Aborting [Undefined neutrino type. Use the --neutrino option]")
unless defined $neutrino;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version; 

$use_valgrind   = 0                            unless defined $use_valgrind;
$arch           = "SL5_64bit"                  unless defined $arch;
$production     = "$genie_version"             unless defined $production;
$cycle          = "01"                         unless defined $cycle;
$batch_system   = "PBS"                     unless defined $batch_system;
$queue          = "prod"                       unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"         unless defined $softw_topdir;
$flux_version   = "10"                         unless defined $flux_version;
$flux_hist_file = "sk_flux_histograms.root"    unless defined $flux_hist_file;
$nevents        = "2000";   
$time_limit     = "05:00:00";
$production_dir = "$softw_topdir/scratch";
$inputs_dir     = "$softw_topdir/data/job_inputs";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$geom_tgt_mix   = "1000080160[0.8879],1000010010[0.1121]";
$xspl_file      = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$flux_file      = "$inputs_dir/t2k_flux/$flux_version/sk/$flux_hist_file";
$job_dir        = "$production_dir/skmc-$production\_$cycle-$neutrino";
$file_prefix    = "genie_sk";

%mcseed_base    = ( 'numu'    => '183221029',
                    'numubar' => '283221029',
                    'nue'     => '383221029',
                    'nuesig'  => '483221029' );
%mcrun_base     = ( 'numu'    => '10000000',
                    'numubar' => '10000000',
                    'nue'     => '10000000',
                    'nuesig'  => '10000000' );
%nu_pdg_code    = ( 'numu'    =>  '14',
                    'numubar' => '-14',
                    'nue'     =>  '12',
                    'nuesig'  =>  '12' );
%flux_hist_name = ( 'numu'    => 'numu_flux',
                    'numubar' => 'numubar_flux',
                    'nue'     => 'nue_flux',
                    'nuesig'  => 'numu_flux' );
#%flux_hist_name = ( 'numu'    => 'h100',
#                    'numubar' => 'h200',
#                    'nue'     => 'h300',
#                    'nuesig'  => 'h100' );

die("** Aborting [Can not find GENIE setup script: ..... $genie_setup]") 
unless -e $genie_setup;
die("** Aborting [Can not find flux file: .............. $flux_file]")   
unless -e $flux_file;
die("** Aborting [Can not find xsec file: .............. $xspl_file]")   
unless -e $xspl_file;

# make the jobs directory
# 
mkpath ($job_dir, {verbose => 1, mode=>0777}); 

# form mcrun, mcseed numbers 
#
$mcrun  = $run + $mcrun_base  {$neutrino}; 
$mcseed = $run + $mcseed_base {$neutrino};

print "@@@ Will submit job with MC run number = $mcrun (seed number = $mcseed)\n";

# get neutrino code and flux histogram name
#
$nu  = $nu_pdg_code    {$neutrino};
$hst = $flux_hist_name {$neutrino};

# form event generation and file conversion commands
#
$fntemplate    = "$job_dir/skjob-$mcrun";
$ghep_file     = "$file_prefix.$production\_$cycle.$neutrino.$mcrun.ghep.root";
$grep_pipe     = "grep -B 50 -A 50 -i \"warn\\|error\\|fatal\"";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_cmd     = "gT2Kevgen -g $geom_tgt_mix -f $flux_file,$nu\[$hst\] -r $mcrun -n $nevents | $grep_pipe &> $fntemplate.evgen.log";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f t2k_tracker -i $ghep_file";

print "@@@ exec: $evgen_cmd \n";

#
# submit
#

# PBS case
if($batch_system eq 'PBS') {
  $batch_script = "$fntemplate.pbs";
  open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
  print PBS "#!/bin/bash \n";
  print PBS "#PBS -N $mcrun\_sk-$production-$cycle \n";
  print PBS "#PBS -l cput=$time_limit \n";
  print PBS "#PBS -o $fntemplate.pbsout.log \n";
  print PBS "#PBS -e $fntemplate.pbserr.log \n";
  print PBS "source $genie_setup \n";
  print PBS "cd $job_dir \n";
  print PBS "export GSPLOAD=$xspl_file \n";
  print PBS "unset GEVGL \n";
  print PBS "export GSEED=$mcseed \n";
  print PBS "$evgen_cmd \n";
  print PBS "$frenm_cmd \n";
  print PBS "$fconv_cmd \n";
  close(PBS);
  `qsub -q $queue $batch_script`;
}

# LSF case
if($batch_system eq 'LSF') {
  $batch_script = "$fntemplate.sh";
  open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
  print LSF "#!/bin/bash \n";
  print LSF "#BSUB-j $mcrun\_sk-$production-$cycle \n";
  print LSF "#BSUB-q $queue \n";
  print LSF "#BSUB-c $time_limit \n";
  print LSF "#BSUB-o $fntemplate.lsfout.log \n";
  print LSF "#BSUB-e $fntemplate.lsferr.log \n";
  print LSF "source $genie_setup \n";
  print LSF "cd $job_dir \n";
  print LSF "export GSPLOAD=$xspl_file \n";
  print LSF "unset GEVGL \n";
  print LSF "export GSEED=$mcseed \n";
  print LSF "$evgen_cmd \n";
  print LSF "$frenm_cmd \n";
  print LSF "$fconv_cmd \n";
  close(LSF);
  `bsub < $batch_script`;
}