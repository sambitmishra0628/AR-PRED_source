#!/usr/bin/perl
# Get the secondary structure assignments for each resdue
# in the PDB structure. Make sure that the dssp path is included
# the PATH variable.

use strict;
use Cwd;
my ($pdb_allatom_file, $feature_file, $outputdir, $logfile) = @ARGV;

eval {
	open LH, ">>", "$logfile" or die "$!\n";
}; if ($@) {
	print "ERROR! Cannot open log file $logfile. $@\n";
	exit;
}	

# Create a separate output subdirectory to store the results from DSSP
my $output_subdir = "$outputdir/dssp_output";
unless (-d $output_subdir) {
	eval {
		mkdir($output_subdir) or die "Cannot create directory $output_subdir $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print LH "ERROR $@\n";
		close LH;
		exit;
	}
}	
my $pdb_allatom_file_fq = "$outputdir/$pdb_allatom_file"; # Get the fully qualified path for the pdb file
my $dssp_output_file = "$output_subdir/dssp.out";

# Execute dssp on the pdb file
print "Running secondary structure calculations for $pdb_allatom_file";
print LH "Running secondary structure calculations for $pdb_allatom_file";
`dssp -i $pdb_allatom_file_fq -o $dssp_output_file`;

# Parse the DSSP output file
my ($ss_seq, $error_flag) = ParseDSSPOutput($dssp_output_file, \*LH);
if ($error_flag) {
	print "ERROR parsing DSSP output!\n";
	print LH "ERROR parsing DSSP output!\n";
	close LH;
	exit;
}

# Add the secondary structure assingment to the feature file
eval {
	open FH, ">", "$outputdir/$feature_file" or die "Cannot create $feature_file $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
my $ss_seq_csv = join(',', split('', $ss_seq)); # convert into comma separated values
print FH "$ss_seq_csv";
close FH;

exit;

#####################################
# ParseDSSPOuput
#	Parse the DSSP output file to get the
#	residue-level secondary structure information.
#
# Input-
# 	dssp_output_file: Name of the output file from DSSP
# 	LH : Handle to the log file
# 
# Output-
# 	parsed secondary structure elements for
# 	each residue
#####################################
sub ParseDSSPOutput  {
	my ($dssp_output_file, $LH) = @_;
	my $errorflag = 0;
	# Map the different secondary structure assignments to
	# numbers.
	my $ss_map = {
		'H' => 1,
		'S' => 2,
		'G' => 3,
		'T' => 4,
		'E' => 5,
		'B' => 6,
		'I' => 7,
		'-' => 8,
	};
	eval {
		open FH, "$dssp_output_file" or die "Cannot open file $dssp_output_file $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print $LH "ERROR $@\n";
		$errorflag = 1;
		return ('', $errorflag);
	}
	my $start_flag = 0;
	my $ss_seq;
	while (my $line = <FH>) {
		next if $line =~ /^\s+$/;
		if ($line =~ /^\s*?#.*?RESIDUE.*?AA.*?STRUCTURE/) { # Marker that our SS info is about to begin
			$start_flag = 1;
			next;
		}
		if ($start_flag == 1) {
			my $sn = substr($line, 0, 5);  $sn =~ s/\s+//g; # serial num
			my $rn = substr($line, 5, 5); $rn =~ s/\s+//g;  # residue num
			my $chain = substr($line, 10, 2); $chain =~ s/\s+//g; # chain
			my $aa = substr($line, 12, 2); $aa =~ s/\s+//g; # amino acid
			my $ss = substr($line, 16, 1); $ss =~ s/\s+//g; # secondary structure assignment
			if ($sn && $rn && $chain && $aa) { # consider the assignment only if these values exist
				if ($ss =~ /\S+/) {
					$ss_seq .= $ss_map->{$ss}; # Get the numeric value for the current secondary structure element
				} else {
					$ss = '-';
					$ss_seq .= $ss_map->{$ss};
				}
			}
		}
	}
	close FH;
	return ($ss_seq, $errorflag);
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

