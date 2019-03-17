#!/usr/bin/perl
###########################################################
# Run calculations for residue conservation scores and write
# residue-level conservation scores into a feature file.
#
# Input arguments:
# 	- pdb c-alpha file name 
# 	- feature file name
# 	- outputdir name
# 	- log file name
#
# Note. The path to blastp, clustalomega and rate4site must
# be included in the PATH variable. This script is intended
# to be invoked from the master_script in the AR_PRED pipeline.
# 
##########################################################

use strict;
use File::Path;
use Cwd;
local $|=1;

### List of pre-defined constants ###
our $max_perc_id = 95; # Maximum percentage identity between the query and subject
our $min_perc_id = 35; # Minimum percentage identity between the query and subject
our $coverage_cutoff = 80; # Minimum coverage we need between the query and subject
our $evalue_cutoff = 0.01; # Maximum e-value for a BLAST hit. Use the default value as by Consurf
our $max_seqs_for_msa = 150; # This is the maximum number of sequences we will use to perform the Multiple Sequence Alignment and then
			     # feed this alignment to rate4site. Since rate4site falls apart when fed with a large MSA file, hence this limit.
our $max_target_seqs = 5000; # The top number of hits from BLAST which we will analyse
# Get the run-time arguments 
my ($pdb_ca_file, $feature_file, $outputdir, $logfile) = @ARGV;

eval {
	open LH, ">>", $logfile or die "$!\n";
}; if ($@) {
	print "ERROR opening log file $logfile $@\n";
	exit;
}
print "\n***Preparing environment for residue conservation calculations***\n";
print LH "\n***Preparing environment for residue conservation calculations***\n";

# Create a separate sub-directory under the outputdir to store
# all the associated files generated during residue conservation
# calculations.
my $output_subdir = "$outputdir/residue_conservation";
unless (-d $output_subdir) {
	eval {
		mkdir ($output_subdir) or die "Cannot create directory $output_subdir $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print LH "ERROR $@\n";
		close LH;
		exit;
	}	
}
# Declare variables and file names
my ($pdbid) = $pdb_ca_file =~ /(^.*?)\.pdb/;
my $blasthitseq = 'blast.hits.seq';
# Create the MSA input and output directory
my $msainputdir = "$output_subdir/ClustalO_input";
eval {
	mkdir($msainputdir) or die "Could not create directory $msainputdir $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
my $msaoutputdir = "$output_subdir/ClustalO_output";
eval {
	mkdir($msaoutputdir) or die " Could not create directory $msaoutputdir $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
# Create an output directory for BLAST results
my $blastoutputdir = "$output_subdir/blast_output";
eval {
	mkdir($blastoutputdir) or die "Could not create directory $msaoutputdir $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
# Create rate4site output directory
my $rate4siteoutputdir = "$output_subdir/rate4site_output";
eval {
	mkdir($rate4siteoutputdir) or die "Could not create directory $rate4siteoutputdir $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}

my ($pdbseq, $errorflag) = GetSeqFrmStruct($pdb_ca_file, $outputdir, \*LH);
if ($errorflag || length($pdbseq) == 0) {
	print "ERROR getting sequence from pdb c-alpha file\n";
	print LH "ERROR getting sequence from pdb c-alpha file\n";
	close (LH);
	exit;
}
#print "Sequence : $pdbseq\n";
my $seqlen = length($pdbseq);
# Write the sequence to a temporary file. We will use this file as input to BLAST.
my $tempseq_file = "$output_subdir/temp.seq";
eval {
	open FH, ">", "$tempseq_file" or die "Could not create file $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
print FH ">$pdbid\n$pdbseq\n";
close FH;

# perform blast
my $blast_output_file = $blastoutputdir . "/" . $pdbid .  "_blast.out";
print "***Running BLAST for $pdb_ca_file***";
print LH "***Running BLAST for $pdb_ca_file***";
my $blastcmd = "blastp -query $tempseq_file -out $blast_output_file -evalue $evalue_cutoff -max_target_seqs $max_target_seqs -db nr -outfmt '7 qseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore'";
my $blastcmd_output = `$blastcmd`;
print LH "BLAST CMD Output $blastcmd_output\n";

# Parse output and Fetch the ids which satisfy the filtering criteria
print "Parsing BLAST output and fetching hit IDs for $pdb_ca_file...";
my ($hitidsfile, $errorflag) = ParseBlast($blastoutputdir,$blast_output_file,$seqlen,$pdbid,\*LH);
if ($errorflag) {
	print "ERROR parsing BLAST results\n";
	print LH "ERROR parsing BLAST results\n";
	close (LH);
	exit;
}
print "Wrote to file $blastoutputdir/$hitidsfile\n";
print LH "Wrote to file $blastoutputdir/$hitidsfile\n";

# Fetch the sequences for the hits
my $hitseqfile = $blastoutputdir . "/" . $pdbid. "_homologs.seq";
print "***Fetching sequences for BLAST Hit IDs***";
print LH "***Fetching sequences for BLAST Hit IDs***";
my $blasthitseq_cmd = "blastdbcmd -entry_batch $hitidsfile -db nr -outfmt %f -target_only -out $hitseqfile";
my $blasthitseqcmd_output = `$blasthitseq_cmd`;
print LH "BLAST DB CMD Output: $blasthitseqcmd_output\n";

# Run CD-HIT on the hit seq file
my $cdhit_seqs = $pdbid . '.cdhitnrseqs.fasta';
print "****Running CD-Hit and obtaining non-redundant sequences for MSA***";
print LH "****Running CD-Hit and obtaining non-redundant sequences for MSA***";
my $cdhit_output_file = "$msainputdir\/$cdhit_seqs";
my $cdhit_cmd = "cd-hit -i $hitseqfile -o $cdhit_output_file -c 0.95"; # use 95% identity cutoff for clustering sequences
my $cdhit_cmd_output = `$cdhit_cmd`;
print LH "CD Hit command output: $cdhit_cmd_output\n";

# Get the final number of homologous sequences 
my $numhomologs = `grep ">" $cdhit_output_file | wc -l`;
$numhomologs =~ s/\s+//g;

# Add the query sequence as well to the beginning of the file
eval {
	open FH, "$cdhit_output_file" or die "Cannot open file $cdhit_output_file $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
my @contents = <FH>;
unshift(@contents,">$pdbid\n$pdbseq\n");
close FH;

# Retain only the number of homologous sequences given by max_seqs_for_msa
my @filtered_contents = FilterMSAInput(\@contents, $max_seqs_for_msa); # Retain only max_seqs_for_msa number of sequences for MSA
if (scalar(@filtered_contents) == 0) {
	print "ERROR No sequences found after filtering MSA\n";
	print LH "ERROR No sequences found after filtering MSA\n";
	close LH;
	exit;
}
my $seqs4msa = $pdbid . '.cdhit_subset_4_msa.fasta';
eval {
	open FH, ">", "$msainputdir\/$seqs4msa" or die "Cannot create file $msainputdir\/$seqs4msa $!\n";
}; if ($@) {
	print "ERROR $@\n";
	print LH "ERROR $@\n";
	close LH;
	exit;
}
print FH @filtered_contents;
close FH;

# Run Clustal Omega on the set of sequences 
my $msa_input_file = "$msainputdir\/$seqs4msa";
my $msa_output_file = $msaoutputdir . "/" . $seqs4msa . ".aln";
print "***Running ClustalO MSA for $seqs4msa***"; 
print LH "***Running ClustalO MSA for $seqs4msa***"; 
my $clustal_cmd = "clustalomega -i $msa_input_file --seqtype=Protein -o $msa_output_file";
my $clustal_output = `$clustal_cmd`;
print LH "Clustal omega output: $clustal_output\n";

# Execute rate4site for the alignment file generated in previous step
print "***Running rate4site for $msa_output_file***\n";
print LH "***Running rate4site for $msa_output_file***\n";

my $rate4siteoutputfile = $rate4siteoutputdir . "/" . $pdbid . ".scores";
# use the same parameters for rate4site as used by Consurf i.e method as Bayesian, 
# evolutionary model as JTT (used as default by rate4site)
my $rate4sitecmd = "rate4site -s $msa_output_file  -ib -o $rate4siteoutputfile";
`$rate4sitecmd`;

# Add conservation scores to feature file
AddConservation2FeatureFile($rate4siteoutputdir, $feature_file, $outputdir, \*LH);
exit;


########################################
# GetSeqFrmStruct
#	Get the amino acid sequence from the pdb C-alpha file
#
# Input Arguments-
# 	- pdb_ca_file: The name of PDB file having coordinates of only 
# 					C-alpha atoms.
# 	- outputdir: Name of the output directory
# 	- LH: Handle to the log file
# Returns-
# 	The single letter amino acid sequence of the protein
######################################## 	

sub GetSeqFrmStruct {
	# Fetch the sequence of the given C_alpha structure 
	my ($pdb_ca_file, $outputdir, $LH) = @_;
	my $pdb_ca_file_fq = "$outputdir/$pdb_ca_file";
	my $errorflag = 0;
	# Define the single letter code for each residue
	my $residue_map = {
		'ALA' => 'A',
		'ARG' => 'R',
		'ASN' => 'N',
		'ASP' => 'D',
		'CYS' => 'C',
		'GLY' => 'G',
		'GLN' => 'Q',
		'GLU' => 'E',
		'HIS' => 'H',
		'ILE' => 'I',
		'LEU' => 'L',
		'LYS' => 'K',
		'MET' => 'M',
		'PRO' => 'P',
		'PHE' => 'F',
		'SER' => 'S',
		'THR' => 'T',
		'TYR' => 'Y',
		'TRP' => 'W',
		'VAL' => 'V',
	};
	eval {
		open FH2, "$pdb_ca_file_fq" or die "Cannot open file $pdb_ca_file_fq $!\n";
	}; if ($@) {
		$errorflag = 1;
		print "ERROR $@\n";
		print $LH "ERROR $@\n";
		return ('', $errorflag);
	}	
	
	my $seq = '';
	print "Fetching sequence for PDB $pdb_ca_file\n";
	print $LH "Fetching sequence for PDB $pdb_ca_file\n";
	while (my $line = <FH2>) {
		next if $line =~ /^\s+$/g; # skip blank lines
		if ($line =~ /^ATOM/) { # Fetch the sequence corresponding only to the given chain
			# We assume that the pdb file is C-alpha only and hence, the residue sequences are in serial NR-order.
			my $res = substr($line,17,3); $res =~ s/\s+//g;
			my $code = $residue_map->{$res};
			$seq .= $code;
		}
	}
	close FH2;
	return $seq;
}

############################################
# ParseBlast
# 	Parse the Blast output and get only those
# 	hit ids which satisfy the criteria for
# 	 - maximum percentage identity
#	 - minimum percentage identity
#	 - coverage
# Input arguments-
# 	blastoutputdir: The directory where all the output files
# 					from BLAST will be placed
# 	blast_output_file: The output file from BLAST having
# 					information on all the homologs collected
# 	seq_len: Length of the query sequence
# 	pdb_id: PDB ID for the structure
# 	LH: The reference to file handle for the log file
#
# Returns-
#	hits_id_file: Name of the file having the IDs of all the
#				  hits collected that satisfy the criteria in place
############################################ 	
sub ParseBlast {
	my ($blastoutputdir, $blast_output_file, $seq_len, $pdb_id, $LH) = @_;
	my $errorflag = 0;
	my $hits_id_file = $blastoutputdir . "/" . $pdb_id . ".ids.txt";
	eval {
		open FH, "$blast_output_file" or die "Cannot read file $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print $LH "ERROR $@\n";
		$errorflag = 1;
		return ('', $errorflag);
	}
	eval {
		open FH2,">","$hits_id_file" or die "Cannot create file $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print $LH "ERROR $@\n";
		$errorflag = 1;
		return ('', $errorflag);
	}
	while (my $line = <FH>) {
		next if $line =~ /^\s+$|^#/; # skip blank lines or comments
		# check for the Gene ID and capture it
		my @line_elements = split(/\s+/, $line);
		my $geneid = $line_elements[1]; 
		my $perc_id = $line_elements[2];
		my $aln_len = $line_elements[3];
			
		# Check the coverage for the alignment. We only consider hits for which
		# the coverage is >= than 80 percent
		my $coverage = ($aln_len/$seq_len)*100;
		
		# Write the IDs which satisfy the cutoff criteria
		if ($perc_id <= $max_perc_id && $perc_id >= $min_perc_id && $coverage >= $coverage_cutoff) {
			print FH2 "$geneid\n";
		}
	}
	close FH;
	close FH2;
	return ($hits_id_file);
}
##################################
# FilterMSAInput
#	Get a subset of the homologs as set by
#	max_seqs_for_msa.
# 
# Input arguments-
# 	contents - array ref for the set of homologs in fasta format
# 	max_seqs_for_msa - Subset of sequences to extract
# 	
# Returns-
# 	Subset of homologs 
################################## 	
sub FilterMSAInput {
	my ($contents, $max_seqs_for_msa) = @_;;
	my @contents = @$contents;
	my @filtered_contents;
	my $numseqs = 0;
	for my $line(@contents) {
		if ($line =~ /^\>/) {
			$numseqs += 1;
			last if $numseqs > ($max_seqs_for_msa+1);
			push (@filtered_contents, $line);
		} else {
			# substitute any non-standard amino acid alphabets other than the
			# 21 amino acids and 'X' by 'X'
			$line =~ s/[^ACDEFGHIKLMNPQRSTVWXY\s]/X/g;
			push (@filtered_contents, $line);
		}
	}
	return @filtered_contents;
}


# Parse the conservation scores generated from rate4site and then add it to the featurefile
sub AddConservation2FeatureFile { 
	my ($scoredir, $featurefile, $outputdir, $LH) = @_;
	my $errorflag = 0;
	eval {
		opendir(DH, $scoredir) or die "Cannot open rate4site scores directory $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print $LH "ERROR $@\n";
		$errorflag = 1;
		return $errorflag;
	}
	my @scorefiles = readdir(DH);
	closedir (DH);

	for my $scorefile(@scorefiles) {
		next if $scorefile !~ /\.scores$/;
		print "Parsing file $scorefile...";
		my ($score_vect, $errorflag) = ParseScoreFile($scoredir,$scorefile,$LH);
		if ($errorflag) {
			print "ERROR parsing score file\n";
			print $LH "ERROR parsing score file\n";
			return $errorflag;
		}
		eval {
			open FH2, ">", "$outputdir/$featurefile" or die "Cannot open file $!\n";
		}; if ($@) {
			print "ERROR $@\n";
			print $LH "ERROR $@\n";
			$errorflag = 1;
			return $errorflag;
		}
		print FH2 "$score_vect\n";
		close FH2;
	}
	return $errorflag;
}

# Parse the Conservation Score generated with Rate4Site
sub ParseScoreFile {
	my ($scoredir, $scorefile, $LH) = @_;
	my $errorflag = 0;
	eval {
		open FH, "$scoredir\/$scorefile" or die "Cannot open file $!\n";
	}; if ($@) {
		print "ERROR $@\n";
		print $LH "ERROR $@\n";
		$errorflag = 1;
		return ([], $errorflag);
	}
	my @score_vect;
	while (my $line = <FH>) {
		next if $line =~ /^\s*?#|^\s+$/; # skip comment lines or blank lines
		$line =~ s/^\s+//g; # remove white space in the beginning of the line
		my @line_elems = split('\s+',$line);
		my $score = $line_elems[2];
		push (@score_vect, $score);
	}
	close FH;
	my $score_vect = join (',',@score_vect);
	return $score_vect;
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
