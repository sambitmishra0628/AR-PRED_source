#!/usr/bin/perl
use strict;

############################################
# This script is intended to be invoked from the master script.
# Given a PDB file, this script will do the following:
# 1. Parse the coordinates from the PDB file for the chain of interest
# 2. Remove anisotropic records
# 3. Remove hetero-atoms
# 4. Remove alternate 'B' locations for residues, if any, and retain only the
#    'A' locations
# 5. Identify any non-standard or modified amino acids that might be present and 
# 	 prompt the user to replace them by their standard names.
# In case of multiple models, this script will consider the first model.
############################################


my ($native_pdb_file, $chain_of_interest, $log_file, $outputdir) = @ARGV;
#our $log_file = 'log.txt'; # store all errors and outputs to this log file
eval {
	open LH, ">>", $log_file or die "Cannot open $log_file $!\n";
}; if ($@) {
	print $@;
	exit;
}	

# Call the individual subroutines. To all the subroutines called, we will
# pass the log file handle to write the progress into the log file. Each subroutine
# will return a boolean variable that will tell us if we encountered any unexpected
# error and terminate execution of the script.
#
# 1. Process pdb and get the records for the chain of interest
print "***Processing pdb file***\n";
print LH "***Processing pdb file***\n";
my ($processed_pdb_array, $errorflag) = process_pdb($outputdir . '/' . $native_pdb_file, \*LH, $chain_of_interest);

# Check if we encountered any error while processing the pdb file
if ($errorflag == 1) {
	print "ERROR! Errors encountered when parsing PDB file! Please check log file!\n";
	print LH "ERROR! Errors encountered when parsing PDB file! Will not proceed with further calculations until the current errors have been fixed!\n";
	close LH;
	exit;
}

# Write the processed pdb coordinates into new file
(my $native_file_prefix) = $native_pdb_file =~ /(^.*?)\.pdb$/;
my $processed_file = $outputdir . "/" . $native_file_prefix . "_processed.pdb";
eval {
	open fh1, ">", "$processed_file" or die "Cannot create file $processed_file $!\n";
}; if ($@) {
	print "ERROR! $@\n";
	print LH "ERROR! $@\n";
	close(LH);
	exit;
} else {	
	print fh1 @$processed_pdb_array;
	close(fh1);
	print "\tWrote processed pdb to $processed_file\n";
	print LH "\tWrote processed pdb to $processed_file\n";
}

# 2. Get the C-alpha coordinates for the PDB and write it into a file
# We will also make a check for the number of residues while writing the 
# C-alpha atoms. If the total number of residues in the PDB file is less the
# required cutoff, then notify the user that the protein is too small.
my $min_res = 50;
print "***Parsing C-alpha coordinates***\n";
print LH "***Parsing C-alpha coordinates***\n"; 
my ($ca_coords_recs, $error_flag) = write_ca_coords($processed_pdb_array, \*LH, $min_res);
if ($error_flag) {
	print "ERROR while parsing C-alpha coordinates! Cannot proceed until current errors are fixed!\n";
	print LH "ERROR while parsing C-alpha coordinates! Cannot proceed until current errors are fixed!\n";
	close(LH);
	exit;
}
my $ca_file = $outputdir . "/" . $native_file_prefix . "_ca.pdb";
eval {
	open fh1, ">", "$ca_file" or die "Cannot create file $ca_file $!\n";
}; if ($@) {
	print "ERROR! $@\n";
	print LH "ERROR! $@\n";
	close(LH);
	exit;
} else {
	print fh1 @$ca_coords_recs;
	close(fh1);
	print "\tWrote C-alpha records to $ca_file\n"; 
	print LH "\tWrote C-alpha records to $ca_file\n";
}


############################################
# Process the pdb file provided by the user by removing anisotropic records,
# retain alternate A locations, remove heteroatoms. If the PDB file has multiple 
# models, then this script will retain the chain of interest only from the first
# model.
# 
# Input arguments -
#	native pdb file
#	log file handle
#	chain of interest
# 
# Returns: 
# 	processed pdb array and error flag
############################################
sub process_pdb {
	my ($native_pdb_file, $LH, $chain_of_interest)  = @_;
	my ($native_file_prefix) = $native_pdb_file =~ /(^.*?)\.pdb$/;
	my @processed_file_contents;
	my $error_flag = 0;
	my @standard_amino_acids = ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLY', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL');
	eval {
		open FH, $native_pdb_file or die "Error! $!";
	}; if ($@) { # check for any possible error while opening the file
		print $LH "$@";
		$error_flag =  1;
		return ([], $errorflag);
	} else {
		my @file_contents = <FH>; # read the file contents into an array. Each line is an element in the array!
		close FH;
		my $default_model = 1;
		my $model_flag = 0;
		my $current_model;
		for my $line_i(@file_contents) {
			# In case the PDB file has multiple models, we will only consider the chain of interest
			# from the first model.
			if ($line_i =~ /^\s*?MODEL(.*?$)/) {
				$current_model = $1;
				$current_model =~ s/\s+//g;
				$model_flag = 1;
				next;
			} elsif ((($model_flag == 1 and $current_model == $default_model) or ($model_flag == 0)) and ($line_i =~ /^ATOM/)) {
				my $atom_name = substr($line_i, 12, 4);
				$atom_name =~ s/\s+//g;
				next if ($atom_name =~ /^H/ && $atom_name ne 'HNCA'); 
				my $alt_loc = substr($line_i,16,1);
				$alt_loc =~ s/\s+//g;
				my $chain = substr($line_i,21,1);
				if (($alt_loc eq '' or $alt_loc eq 'A') and ($chain eq $chain_of_interest)) { # only consider those records not having any alternate location indicators 
																							 # or 'A' locators for those with alternate locations and belonging to the chain
																							 # of interest.
					# First make sure that the residue name conforms to the standard set of
					# 20 amino acids. Otherwise, instruct the user to change the residue name in the
					# PDB file itself.
					my $res_name = substr($line_i, 17,3);
					if (grep {$_ eq $res_name} @standard_amino_acids) {
						push (@processed_file_contents, $line_i);
					} else {
						print "Non-standard amino acid found in PDB file!\n$line_i Please make sure all residues are standard amino acids in the PDB file!\n";
						print $LH "Non-standard amino acid found in PDB file!\n$line_i Please make sure all residues are standard amino acids in the PDB file!\n";
						$error_flag=1;
						return ([], $error_flag);
					}
				} else {
					next;
				}

			} elsif ($line_i =~ /^\s*?ENDMDL/ and $current_model == $default_model) {# if we reach the end of the first model let's quit parsing
				print $LH "PDB file contains multiple models. Proceeding with first model!\n";
				last;
			}
		}
	}
	if (scalar(@processed_file_contents) <= 1) { # In case we are not able to parse the records for the given chain
		print $LH "Error! parsing PDB contents! No ATOM records found for $chain_of_interest\n";
		$error_flag = 1;
	}	
	return (\@processed_file_contents, $error_flag);
}
#######################################################
# Get the C-alpha coordinates for the PDB chain of interest
# Input arguments:
# 	processed_pdb_array: Array of all processed all-atom records
# 	LH: log file handle
# 	min_res: Minimum number of residues required.
#
# Returns:
# 	The records for C-alpha coordinates and error flag
#######################################################

sub write_ca_coords {
	my ($processsed_pdb_array, $LH, $min_res) = @_;
	my $error_flag = 0;
	my @ca_records = ();
	for my $rec(@$processsed_pdb_array) {
		my $atom_name = substr($rec, 12,4);
		$atom_name =~ s/\s+//g;
		if ($atom_name eq 'CA') {
			push @ca_records, $rec;
		}	
	}
	
	if (scalar(@ca_records)<$min_res) {
		print "ERROR! Protein too short. Must have atleast $min_res residues with C-alpha atoms for calculations.\n";
		print $LH "ERROR! Protein too short. Must have atleast $min_res residues with C-alpha atoms for calculations.\n";
		$error_flag = 1;
	}
	return (\@ca_records, $error_flag); 
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

