#!/usr/bin/perl
use strict;
use Cwd;

# The master script!!! It's a parent script which will be used to invoke
# all other processing/calculation scripts. Make sure all the required environment
# variables are set prior to running this script. 

my $USAGE = <<"USAGE";
Usage: perl $0 pdb_file chain outputdir pred_type activesite_file

pdb_file: Name of the pdb file if file is present in current directory otherwise,
absolute or relative path must be included as part of filename (Required)

chain: The chain for which predictions are to be made (Required)

outputdir: The outputdir where the prediction output will be stored (Required)

pred_type: Set to 1 when predicting active site residues,
		 and 2 when predicting allosteric site residues. (Required)

activesite_file: CSV file having the list of active site residues. Required iff
		pred_type is 2
			
USAGE

my ($pdb_file, $chain, $outputdir, $pred_type, $activesite_file) = @ARGV;

die "$USAGE\n" unless ($pdb_file && $chain && $outputdir && $pred_type);

if ($pred_type != 1 && $pred_type != 2) {
	die "ERROR with argument pred_type\n$USAGE";
} elsif ($pred_type == 2 && not $activesite_file) {
	die "ERROR! Active site file is required with pred_type 2\n$USAGE";	
}

# Make sure all the environment variables are set. List of environment variables to be set prior to running the
# script are as follows.
# 1. ARPRED_HOME
# 2. MATLAB_PATH
unless ($ENV{ARPRED_HOME} && $ENV{MATLAB_PATH} && $ENV{BLASTDB} && $ENV{NACESS_PATH}) {
	die "The environment variables ARPRED_HOME,MATLAB_PATH,BLASTDB and NACESS_PATH must be set before running this code. Please refer to the instructions in \'Running AR-PRED.txt\'!";
}

$outputdir =~ s/\/$//g; # remove any forward slashes that user might have accidentally given

unless (-d $outputdir) {
	`mkdir $outputdir`;
}
# chdir into outputdir to get the full path of the outputdir
my $current_dir = getcwd();
chdir ($outputdir);
my $outputdir_fq = getcwd();
chdir ($current_dir);

my $log_file = "$outputdir\/log.txt";
open (FH,">", $log_file) or die "Cannot create log file $!\n";
close(FH);

# 1. Process the PDB file
die "ERROR! PDB file needs to end with extension .pdb\n" if $pdb_file !~ /\.pdb$/;

my $pdb_file_prefix = $pdb_file;
$pdb_file_prefix =~ s/\.pdb$//;

my $pdb_allatom_file = $pdb_file_prefix . "_processed.pdb"; # This is the naming convention used by verify_and_process_pdb 
my $pdb_ca_file = $pdb_file_prefix . "_ca.pdb";
if (-e "$outputdir_fq\/$pdb_allatom_file" && -e "$outputdir_fq\/$pdb_ca_file") {
	print "Skipping processing of pdb file as processed files are already present!\n";
} else {
	print "Processing PDB...\n"; 
	my $script_name = "$ENV{ARPRED_HOME}/pdb_file_processing_scripts/verify_and_process_pdb.pl"; 
	my $output = `perl $script_name $pdb_file $chain $log_file $outputdir_fq`;

	if ($output =~ /ERROR/) {
		die "Error encountered in $script_name! Please see log file for details.";
	} else {
		print "Done!\n";
	}
}

# 2. Run Solvent accessibility calculations
my $sasa_feature_file = 'sasa_features.csv';
if (-e "$outputdir_fq\/$sasa_feature_file") {
	print "Skipping calculations for solvent accessibility as feature file is already present\n";
} else {
	print "Running calculations for Solvent Accessibility...\n";
	my $script_name = "$ENV{ARPRED_HOME}/solvent_accessibility/calc_sasa.pl";
	my $output = `perl $script_name $pdb_allatom_file $sasa_feature_file $outputdir_fq $log_file`;
	
	if ($output =~ /ERROR/) {
		die "Error encountered in $script_name! Please see log file for details.";
	}
	else {
		print "Done!\n";
	}	
}
# 3. Run FPocket calculations
my $fpocket_feature_file = 'pocket_features.csv';
if (-e "$outputdir_fq\/$fpocket_feature_file") {
	print "Skipping calculations for pocket features as feature file is already present\n";
} else {	
	print "Running Pocket calculations...\n";
	my $script_name = "$ENV{ARPRED_HOME}/pocket_residues/GetPocketResidues.pl";
	my $output = `perl $script_name $pdb_allatom_file $fpocket_feature_file $outputdir_fq $log_file`;

	if ($output =~ /ERROR/) {
		die "Error encountered in $script_name! Please see log file for details.";
	}
	else {
		print "Done!\n";
	}	
}

# 4. Run Residue Conservation Calculations
my $conservation_feature_file = 'conservation.csv';
if (-e "$outputdir_fq\/$conservation_feature_file") {
	print "Skipping calculations for conservation scores as feature file is already present\n";
} else {	
	my $script_name = "$ENV{ARPRED_HOME}/residueconservation/residueconservation.pl";
	print "Running calculations for residue conservation...\n";
	my $output = `perl $script_name $pdb_ca_file $conservation_feature_file $outputdir_fq $log_file`;
	if ($output =~ /ERROR/) {
		die "Error encountered in $script_name! Please see log file for details.";
	}
	else {
		print "Done!\n";
	}	
}

# 5. Run DSSP Secondary Structure Calculations
my $dssp_feature_file = 'sse.csv';
if (-e "$outputdir_fq\/$dssp_feature_file") {
	print "Skipping calculations for secondary structure elements as feature file is already present\n";
} else { 	
	my $script_name = "$ENV{ARPRED_HOME}/secondary-structure/get_sec_struc_elems.pl";
	print "Running calculations for secondary structure elements with DSSP...\n";
	my $output = `perl $script_name $pdb_allatom_file $dssp_feature_file $outputdir_fq $log_file`;
	if ($output =~ /ERROR/) {
		die "Error encountered in $script_name! Please see log file for details.";
	}
	else {
		print "Done!\n";
	}	
}

# 6. Run calculations for protein dynamics features with matlab
my $matlab = $ENV{'MATLAB_PATH'}; # The MATLAB_PATH environmental variable needs to be set to the matlab executable

my $dl = '\''; # The string arguments passed to a matlab command interface must be enclosed by a single quote
my $cmd_sec1 = 'path_arpred = getenv(\'ARPRED_HOME\'); all_paths = genpath(path_arpred); addpath(all_paths);';
my $cmd_sec2 = 'matlab_feature_calc_wrapper(' . $dl . $pdb_ca_file . $dl . ', ' . $dl . $outputdir_fq . $dl . ',' . $pred_type .  ',' . ($pred_type == 2 ? $dl . $activesite_file . $dl : $dl.$dl) . '); exit;';
my $matlab_cmd = "$matlab -nodisplay -nojvm -r " . '"' . $cmd_sec1 . $cmd_sec2 . '"';
print "Running calculations for protein dynamics features...\n";
my $output = `$matlab_cmd`;
if ($output =~ /ERROR/) {
	die "Error encountered in protein dynamics feature calculations! Please see log file for details.";
} else {
	print "Done!\n";
}

# 7. Run predictions for the given prediction type
if ($pred_type == 1) {
	my $cmd_sec1 = 'path_arpred = getenv(\'ARPRED_HOME\'); all_paths = genpath(path_arpred); addpath(all_paths);';
	my $func = 'PredictActiveSiteResidues(' . $dl . $outputdir_fq . $dl . ',' . $dl . $pdb_ca_file . $dl . ');exit;';
	my $matlab_cmd = "$matlab -nodisplay -nojvm -r " . '"' . $cmd_sec1 . $func . '"';
	print "Running predictions for Active site residues...\n";
	my $output = `$matlab_cmd`;
	if ($output =~ /ERROR/) {
		die "Error encountered in active site predictions! Please see log file for details.";
	} else {
		print "Done!\n";
		print "Wrote predictions to $outputdir_fq/activesite_predictions.csv\n";
	}
} else {
	my $cmd_sec1 = 'path_arpred = getenv(\'ARPRED_HOME\'); all_paths = genpath(path_arpred); addpath(all_paths);';
	my $func = 'PredictAllostericSiteResidues(' . $dl . $outputdir_fq . $dl . ',' . $dl . $pdb_ca_file . $dl . ');exit;';
	my $matlab_cmd = "$matlab -nodisplay -nojvm -r " . '"' . $cmd_sec1 . $func . '"';
	print "Running predictions for Allosteric site residues\n";
	my $output = `$matlab_cmd`;
	if ($output =~ /ERROR/) {
		die "Error encountered in allosteric site predictions! Please see log file for details.";
	} else {
		print "Done!\n";
		print "Wrote predictions to $outputdir_fq/allostericsite_predictions.csv\n";
	}
}
exit;

### LGPL
#    This file is part of AR-PRED: Active and Regulatory site Prediction
#
#    AR-PRED is a free software. You can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AR-PRED is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    and the GNU Lesser General Public License along with the AR-PRED source code.
#    If not, see <http://www.gnu.org/licenses/>.
