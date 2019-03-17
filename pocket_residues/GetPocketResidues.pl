#!/usr/bin/perl
############################################
# Identify residues which are in pockets (boolean feature : 0 or 1)
# with FPocket for the given PDB all-atom file and then write this
# as a feature into the feature file.
# Input arguments-
# 	pdb file
# 	feature file
# 	outputdir
# 	log file
############################################

use strict;
use Cwd;
local $|=1; # Flush buffer instantly

my ($pdb_allatom_file, $featurefile, $outputdir, $logfile) = @ARGV;
my $currentdir = cwd();
my $featurefile = $outputdir . '/' . $featurefile;

eval {
	open LH, ">>", "$logfile" or die "$!\n";
}; if ($@) {
	print "ERROR! Cannot open log file $logfile. $@\n";
	exit;
}	
# create a temporary directory to copy the PDB file and store the results from Fpocket
my $tempdir = "$outputdir\/fpocket_results/";
if (-d "$tempdir") {
	`rm -rf $tempdir`;
}
mkdir($tempdir) or die "Cannot create directory $!\n";
my ($pdb) = $pdb_allatom_file =~ /(^.*?)\.pdb/;
my $destfile = $tempdir . '/' . $pdb_allatom_file;
`cp "$outputdir\/$pdb_allatom_file" $destfile`;

# Run FPocket algorithm
print "\n***Running Fpocket for file $pdb_allatom_file***\n";
print LH "\n***Running Fpocket for file $pdb_allatom_file***\n";

my $fpocket_output = `fpocket -f $destfile`;
print LH "\tFpocket Output: $fpocket_output\n";

# Parse the pocket output directory from FPocket. FPocket creates an output directory 
# with the naming convention pdbid_out/pocket. 
my $fpocket_output_dir = $tempdir . "/" . $pdb . "_" . 'out' . '/pockets'; # The directory where we have the pocket specific PDB files
print "\tParsing results from FPocket\n";
print LH "\tParsing results from FPocket\n";
my ($pocketres) = IdentifyPocketResidues($fpocket_output_dir);
my @pocketres = @$pocketres;

my ($pdbindices) = GetPDBIndex("$outputdir\/$pdb_allatom_file");
my @pdbindices = @$pdbindices;

my @feature_vect;
for my $res (@pdbindices) {
	if (grep {$_ == $res} @pocketres) {
		push (@feature_vect, 1);
	} else {
		push (@feature_vect, 0);
	}		
}	

# Write the feature into the feature file
eval {
	open FH, ">", "$featurefile" or die "Cannot open file $featurefile $!\n";
}; if ($@) {
	print "ERROR creating Fpocket feature file $@\n";
	print FH "ERROR creating Fpocket feature file $@\n";
	exit;
}
print FH join(",",@feature_vect),"\n";
close FH;
print "\tWrote Pocket Features to $featurefile!\n";
print LH "\tWrote Pocket Features to $featurefile!\n";

# Get the PDB Index from the all-atom pdb file
sub GetPDBIndex {
	my ($sourcefile) = @_;
	my @pdbindices;
	open FH1, "$sourcefile" or die "ERROR Cannot open file $!\n";
	while (my $line = <FH1>) {
		if ($line =~ /^ATOM/) {
			my $resn = substr($line, 22,4);
			$resn =~ s/\s*//g;
			unless (grep {$_ == $resn} @pdbindices) {
				# include the index if not already included
				push (@pdbindices, $resn);
			}	
		}
	}
	close FH1;
	return (\@pdbindices);
}

# Parse the FPocket results to identify residues which lie in pockets
sub IdentifyPocketResidues{
	my ($dir) = @_;
	opendir(DH, $dir) or die "ERROR Cannot open directory $!\n";
	my @files = readdir(DH);
	closedir(DH);
	my @pocket_res; # Get the set of residues which are involved in constituting atleast one pocket
	for my $file (@files) {
		next if $file !~ /\.pdb$/;
		my $pdbsourcefile = $dir . '/' . $file;
		my $pocketscore;
		my $pocket_ind;
		open FH, "$pdbsourcefile" or die "Cannot open file $!\n";
		while (my $line = <FH>) {
			if ($line =~ /^ATOM/) {
				my $resn = substr($line,22,4);
				$resn =~ s/\s+//g;
				push (@pocket_res, $resn) unless grep {$_ == $resn} @pocket_res; # Get the unique residue PDB index included in pocket
			}									
		}
		close FH;
	}
	return \@pocket_res;
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

