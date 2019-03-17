#!/usr/bin/perl
###########################################
# Run calculations for residue-level solvent accessibility 
# for the given pdb file using NACCESS.
#
# Input arguments:
# 	- pdb all-atom file name (as in the
# 		output directory)
# 	- feature file name
# 	- outputdir name
# 	- log file name
# The NACCESS path must be included in the PATH variable	
###########################################

use strict;
use Cwd;
local $|=1; # flush buffer

my ($pdb_allatom_file, $feature_file, $outputdir, $logfile) = @ARGV;
eval {
	open LH, ">>", "$logfile" or die "$!\n";
}; if ($@) {
	print "ERROR! Cannot open log file $logfile. $@\n";
	exit;
}	
my ($pdbid) = $pdb_allatom_file =~ /(^.*?)\.pdb/; 
my $naccess_path = $ENV{NACCESS_PATH};

# Copy the pdb file from the current directory to the naccess path and then cd to the path. 
# This step is sometimes needed to make NACCESS run.
`cp "$outputdir\/$pdb_allatom_file" $naccess_path`;
my $current_dir = getcwd();
chdir ($naccess_path);

print "\n***Running Accessibility Calculations for $pdb_allatom_file***\n";
print LH "\n***Running Accessibility Calculations for $pdb_allatom_file***\n";
my $naccess_output = `naccess $pdb_allatom_file`; # run naccess
print LH "\tNACCESS Message: $naccess_output\n";
if ($naccess_output !~/Running Accessibility calculations/) { # Means NACCESS has not run correctly
	print "ERROR executing NACCESS. Please check log file\n";
	print LH "ERROR executing NACCESS. Please check log file\n";
	exit;
}
eval {
	unlink ($pdb_allatom_file) or die "Cannot remove file $!\n";
}; if ($@) {
	print "Warning! Unable to delete copied $pdb_allatom_file from naccess directory! $@\n";
	print LH "Warning! Unable to delete copied $pdb_allatom_file from naccess directory! $@\n";
}	

# Naccess creates two files: relative solvent accessibility (rsa) and absolute 
# solvent accessibility (asa). Move the files to the output directory and then
# return to the working dir.
`mv "$pdbid.rsa" "$outputdir\/"`;
`mv "$pdbid.asa" "$outputdir\/"`;
`mv "$pdbid.log" "$outputdir\/${pdbid}_naccess.log"`;
chdir ($current_dir); # return to the current working directory

# Parse the .rsa file and write the extracted features into the feature file
my $rsafile = $pdbid . ".rsa";
my $error_flag = AddAccessibilityFeature($outputdir, $rsafile, $feature_file, \*LH);
if ($error_flag) {
	print "ERROR adding parsed NACCESS features to feature file! Please check log file for details\n";
	print LH "ERROR adding parsed NACCESS features to feature file! Cannot proceed with further calculations\n";
} else {
	print "\tWrote extracted SASA features to $feature_file\n";
	print LH "\tWrote extracted SASA features to $feature_file\n";
}
exit;



# Add the calculated solvent accessibilities from naccess to the feature file 
sub AddAccessibilityFeature {
	my ($outputdir, $rsafile, $featurefile, $LH) = @_;
	
	# variables to store the relative and absolute accessibilities
	my @asa_rel_allatom;
	my @asa_rel_sidechain;
	my @asa_rel_mainchain;
	my @asa_rel_nonpolar;
	my @asa_rel_polar;
	
	my @asa_abs_allatom;
	my @asa_abs_sidechain;
	my @asa_abs_mainchain;
	my @asa_abs_nonpolar;
	my @asa_abs_polar;

	my $error_flag = 0; # Initialize the error flag to 0
	# Parse accessibility file
	eval {
		open FH, "$outputdir/$rsafile" or die "Cannot open file $!\n";
	}; if ($@) {
		print "ERROR opening $rsafile!\n";
		print $LH "ERROR opening $rsafile!\n";
		$error_flag = 1;
		return $error_flag;
	}
	while (my $line = <FH>) {
		next if $line !~ /^RES/;
		chomp $line; 
		my $line_len = length($line);
		my $line_substr = substr($line,13,$line_len-13); # let's consider only the last 10 cols corresponding to the ASA values
		if ($line_substr =~ /^\s+/) {
			$line_substr =~ s/^\s+//g;
		}
		my @line_elems = split(/\s+/, $line_substr);
		push (@asa_abs_allatom, $line_elems[0]);
		push (@asa_rel_allatom, $line_elems[1]);
		push (@asa_abs_sidechain, $line_elems[2]);
		push (@asa_rel_sidechain, $line_elems[3]);
		push (@asa_abs_mainchain, $line_elems[4]);
		push (@asa_rel_mainchain, $line_elems[5]);
		push (@asa_abs_nonpolar,  $line_elems[6]);
		push (@asa_rel_nonpolar,  $line_elems[7]);
		push (@asa_abs_polar, $line_elems[8]);
		push (@asa_rel_polar, $line_elems[9]);
	}
	close FH;
	
	# write it into the feature file
	eval {
		open FH1, ">", "$outputdir\/$featurefile" or die "Could not open file in append mode $!\n";
	}; if ($@) {
		print "ERROR opening $featurefile $@\n";
		print $LH "ERROR opening $featurefile $@\n";
		$error_flag = 1;
		return $error_flag;
	}
	print FH1 join(",", @asa_abs_allatom),"\n";
	print FH1 join(",", @asa_rel_allatom),"\n";
	print FH1 join(",", @asa_abs_sidechain),"\n";
	print FH1 join(",", @asa_rel_sidechain),"\n";
	print FH1 join(",", @asa_abs_mainchain),"\n";
	print FH1 join(",", @asa_rel_mainchain),"\n";
	print FH1 join(",", @asa_abs_nonpolar),"\n";
	print FH1 join(",", @asa_rel_nonpolar),"\n";
	print FH1 join(",", @asa_abs_polar),"\n";
	print FH1 join(",", @asa_rel_polar),"\n";
	close FH1;
	return $error_flag;
}

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
