#!/usr/bin/env perl

# Runs the software rRNA_predition to find all the rRNA segments in a metagenome
# sample set. Need to have qiime/1.8.0 loaded and seqtk

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use File::Basename;
use BioUtils::FastaIO;

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($exe, $sample_dir, $split_num, $file_ident, $seqtk, $help, $man);

my $options_okay = GetOptions (
    "exe|e:s" => \$exe,
    "sample_dir|d:s" => \$sample_dir,
	"num|n:i" => \$split_num,
    "file_ident|fi:s" => \$file_ident,
	"seqtk:s" => \$seqtk,
    "help|h" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# set up the logging environment
my $logger = get_logger();

# check for input errors
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();

########
# MAIN #
########

# Check the file
my $SAM;
opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");
    
    #go into each dir and get the file
    opendir my $DIR, "$sample_dir/$dir" or
		$logger->logide("Cannot open dir: $sample_dir");
	
	while ( my $file = readdir($DIR) ) {
		if ( $file =~ m/^\./ ) {next;} # skip hidden files
        if ( $file !~ qr/$file_ident/i ) {next;}
        $logger->debug( "Looking at fasta: $sample_dir/$dir/$file" );
		
		#Check if file needs to be unzipped or changed to fasta
		my ( $prefix, $new_file ) = create_fasta_get_prefix( $file, "$sample_dir/$dir", $seqtk );
		
		# generate the split command and run
		my $split_command = "split -l $split_num ";
		$split_command .= "$sample_dir/$dir/$new_file $sample_dir/$dir/$prefix\.";

		# running split command
		$logger->info("Running split command: $split_command");
		`$split_command`;
	}
	closedir($DIR);
}
closedir($SAM);

# create the directories and move the matching file into its directory
opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");
	
	# get the prefix of each fasta in $dir
	opendir my $DIR, "$sample_dir/$dir" or
		$logger->logide("Cannot open dir: $sample_dir");
	
	while ( my $file = readdir($DIR) ) {
		if ( $file =~ m/^\./ ) {next;} # skip hidden files
		if ( $file =~ m/fast(a|q)/ ) {next;} # skip fasta files
		if ( -d "$sample_dir/$dir/$file" ) {next;}
		$logger->debug("Looking at split fasta: $sample_dir/$dir/$file");
		
		my $mv_command = "";
		
		# mv the file to tmp
		$mv_command .= "mv $sample_dir/$dir/$file $sample_dir/$dir/tmp; ";
		
		# make a directory
		$mv_command .= "mkdir $sample_dir/$dir/$file; ";
		
		# move tmp back into the dir
		$mv_command .= "mv $sample_dir/$dir/tmp $sample_dir/$dir/$file/$file; ";
		
		# running all the mv commands
		$logger->info("Running mv command: $mv_command");
		`$mv_command`;
	}
	closedir($DIR);	
}
closedir($SAM);

# start a rRNA_pred job for each split directory
opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");
	
	# get the prefix of each fasta in $dir
	opendir my $DIR, "$sample_dir/$dir" or
		$logger->logide("Cannot open dir: $sample_dir");
	
	while ( my $split_dir = readdir($DIR) ) {
		if ( $split_dir =~ m/^\./ ) {next;} # skip hidden dirs
		if ( ! -d "$sample_dir/$dir/$split_dir" ) {next;} # skip non-dirs
		if ( $split_dir !~ qr/$file_ident/i ) {next;}
		$logger->debug("Looking at split fasta: $sample_dir/$dir/$split_dir");
		
		# build the rRNA_prediction command
		my $pre_command = "bsub -q day -M 8 -J pred ";
		$pre_command .= "-o $sample_dir/$dir/$split_dir/lsf.out ";
		$pre_command .= "-e $sample_dir/$dir/$split_dir/lsf.err ";
		$pre_command .= "$exe ";
		$pre_command .= "-i $sample_dir/$dir/$split_dir/ ";
		$pre_command .= "-o $sample_dir/$dir/$split_dir/ ";
		
		# run the rRNA_prediction command
		$logger->info("Running prediction command: $pre_command");
		`$pre_command`;
	}
	closedir($DIR);	
}
closedir($SAM);

# stall until all the rRNA_pred jobs are finished
while ( stall("pred", $sample_dir) ) {
	$logger->info("Stalling");
	sleep 100;
}

# cat all the rRNA coord files together
opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");
	
	# set up the cat command
	my $cat_command = "cat $sample_dir/$dir/*/*.coord > $sample_dir/$dir/rRNA_pred.coord";
	my $seq_command = "cat $sample_dir/$dir/*/*.seq > $sample_dir/$dir/rRNA_pred.seq";
	
	# run the cat command
	$logger->info("Running cat commands: $cat_command");
	`$cat_command`;
	`$seq_command`;
}

# make the sample count tables and final count table
my $final_count_tbl = "$sample_dir/all_rRNA_pred.count";
open my $OUT, ">", $final_count_tbl or
	$logger->logdie("Cannot open final count tbl file: $final_count_tbl");

my @headers = ("sample", "kingdom", "seq_type", "count");
print $OUT join("\t", @headers), "\n";

opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");
	
	# make count command
	my $coord_file = "$sample_dir/$dir/rRNA_pred.coord";
	my $count_file = "$sample_dir/$dir/rRNA_pred.count";
	my $count_command = "awk \'{print \$8}\' $coord_file | grep -v \"gene\" | sort | uniq -c > $count_file";
	
	# run the count command
	$logger->info("Running count command: $count_command");
	`$count_command`;
	
	# store the count data
	open my $COU, "<", $count_file or
		$logger->logdie("Cannot open count file: $count_file");
	
	foreach my $line ( <$COU> ) {
		chomp $line;
		
		if ( $line =~ m/\s*(\d+)\s+(\w+):(\w+)/ ) {
			print $OUT $dir, "\t", $2, "\t", $3, "\t", $1, "\n";
		}
		elsif ( $line =~ m/\s*(\d+)\s+(Eukaryotic)(18S_rRNA)/ ) {
			# there is a formating error made my rRNA_pred software where
			# they forget to include the ":" between Eukaroytic and 18S_rRNA.
			# this if statement handles that case.
			print $OUT $dir, "\t", $2, "\t", $3, "\t", $1, "\n";
		}
	}
	
	close($COU);
}

close($OUT);

#Run assign taxonomy on output
### assign taxonomy
#  1. split the rRNA_pred output seq file into bacterial 16 and eukaryotic 28S
#  2. run assign_taxonomy.py script in qiime
$logger->debug("Start assign_taxonomy block");

# open the input file
opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");

	my $seq_file = "$sample_dir/$dir/rRNA_pred.seq";
	my $in_fasta = BioUtils::FastaIO->new({stream_type => '<', file => $seq_file});
	
	# open the output files
	my $bac_file = "$sample_dir/$dir/rRNA_pred.seq.16S";
	my $euk_file = "$sample_dir/$dir/rRNA_pred.seq.28S";
	my $bac_out = BioUtils::FastaIO->new({stream_type => '>', file => $bac_file});
	my $euk_out = BioUtils::FastaIO->new({stream_type => '>', file => $euk_file});
	
	# go through each sequence in the input file
	# and seperate out the 16S and 28S sequences
	while ( my $seq = $in_fasta->get_next_seq() ) { 
		my $header = $seq->get_header();
		if ( $header =~ m/Bacterial:16S_rRNA/ ) { 
			$bac_out->write_seq($seq);
		}   
		elsif ( $header =~ m/Eukaryotic:28S_rRNA/ ) { 
			$euk_out->write_seq($seq);
		}   
	}
	
	# check that the output files were made
	if ( ! -s "$sample_dir/$dir/rRNA_pred.seq.16S" ) {
		$logger->warn("Cannot find rRNA_pred_16S_bacterial.seq file for $dir");
		next;
	}
	if ( ! -s "$sample_dir/$dir/rRNA_pred.seq.28S" ) {
		$logger->warn("Cannot find rRNA_pred_28S_eukaryote.seq file for $dir");
		next;
	}
	
	# Submit assign_taxonomy.py job
	my $command = "";
	$command = "bsub -q week -J assign_tax -o $sample_dir/$dir/lsf.out -e $sample_dir/$dir/lsf.err ";
	$command .= "assign_taxonomy.py ";
	$command .= "-i $sample_dir/$dir/rRNA_pred.seq.16S ";
	$command .= "-o $sample_dir/$dir/RDPClassifier_16S/ ";
	$command .= "-r /nas02/home/y/o/yourston/perl5/lib/perl5/MT_MTToolbox/data/97_otus.fasta ";
	$command .= "-t /nas02/home/y/o/yourston/perl5/lib/perl5/MT_MTToolbox/data/97_otu_taxonomy.tab ";
	$logger->info("Submitting 16S assign_taxonomy command: $command");
	`$command`;
	
	$command = "bsub -M 30 -q week -J assign_tax -o $sample_dir/$dir/lsf.out -e $sample_dir/$dir/lsf.err ";
	$command .= "assign_taxonomy.py ";
	$command .= "-i $sample_dir/$dir/rRNA_pred.seq.28S ";
	$command .= "-o $sample_dir/$dir/RDPClassifier_28S/ ";
	$command .= "-r /netscr/yourston/kingdom_proportions/silva_lsu_ref/SILVA_123.1_LSURef_tax_silva_trunc_formated.fasta ";
	$command .= "-t /netscr/yourston/kingdom_proportions/silva_lsu_ref/SILVA_123.1_LSURef_tax_silva_trunc_formated.fasta.tax ";
	$logger->info("Submitting 28S assign_taxonomy command: $command");
	`$command`;
}

while ( stall("assign_tax", $sample_dir) ) {
	$logger->info("Assign_Tax_Stalling");
	sleep 100;
}

# make the final taxonomy table where all samples are in the same file
# open the input file
opendir $SAM, "$sample_dir" or
	$logger->logdie("Cannot open --sample_dir ($sample_dir)");
	
my $tax_16S_file = "$sample_dir/16S_taxonomy_tab.txt";
open my $TAX16, ">", $tax_16S_file
	or $logger->logdie("Cannot open taxonomy_tab.txt file: $tax_16S_file");
print $TAX16 "sample\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";

my $tax_28S_file = "$sample_dir/28S_taxonomy_tab.txt";
open my $TAX28, ">", $tax_28S_file
	or $logger->logdie("Cannot open taxonomy_tab.txt file: $tax_28S_file");
print $TAX28 "sample\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";
	
while ( my $dir = readdir($SAM) ) {
	if ( $dir =~ m/^\./ ) {next;} # skip hidden files
	if ( ! -d "$sample_dir/$dir" ) {next;} # skip non-directories
	$logger->debug("Looking at dir: $dir");
	
	my $file = "$sample_dir/$dir/RDPClassifier_16S/rRNA_pred.seq_tax_assignments.txt";
	open my $IN, "<", $file
		or $logger->logdie("Cannot open taxonomy file: $file");
	
	foreach my $line ( <$IN> ) {
		chomp $line;
		
		my @vals = split(/\t/, $line);
		my @tax_vals = split(/; /, $vals[1]);
		unshift(@tax_vals, $dir);
		if ( $tax_vals[1] eq "Unassigned" ) {
			$tax_vals[1] = "Unassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned";
		}
		my $out_line = join("\t", @tax_vals) . "\n";
		$out_line =~ s/-/_/g;
		print $TAX16 $out_line;
		
	}
	close($IN);
	
	$file = "$sample_dir/$dir/RDPClassifier_28S/rRNA_pred.seq_tax_assignments.txt";
	open $IN, "<", $file
		or $logger->logdie("Cannot open taxonomy file: $file");
	
	foreach my $line ( <$IN> ) {
		chomp $line;
		
		my @vals = split(/\t/, $line);
		my @tax_vals = split(/; /, $vals[1]);
		unshift(@tax_vals, $dir);
		print $TAX28 join("\t", @tax_vals), "\n";
		
	}
	close($IN);
}



        
########
# Subs #
########
sub create_fasta_get_prefix {
	my ( $file_name, $path, $seqtk ) = @_;
	my $name;
	my $new_name;
	if ( $file_name =~ qr/\.gz/i ) {
		my $cmd = "gunzip $path/$file_name";
		system( $cmd );
		$file_name = fileparse( $file_name, qr/\.gz/ );
	}
	if ( $file_name !~ qr/\.fasta/i ) {
		if ( $file_name =~ qr/\.fastq/i ) {
			$logger->debug("Creating fastq file into fasta file");
			
			$name = fileparse( $file_name, qr/\.fastq/ );
			$new_name = "$name.fasta";
			my $cmd = "bsub -J fq -o $path/lsf.out ";
			$cmd .= "'$seqtk seq -A $path/$file_name > $path/$new_name'";
			system ( $cmd );
			while ( stall("fq", $sample_dir) ) {
				$logger->info("Stalling for fq->fa");
				sleep 100;
			}
			if ( -z $path/$new_name ) {
				$logger->logdie("Fastq->Fasta conversion did not work. Make sure path to seqtk is given and correct");
			}
		}
		else {
			carp("$file_name has an unknown file extension. Must be fasta or fastq");
		}
		return ( $name, $new_name );
	}
	else {
		$name = fileparse( $file_name, qr/\.fastq/ );
		return ( $name, $file_name );
	}
}

sub check_params {
	# check for required variables
	if ( ! defined $exe) { 
		pod2usage(-message => "ERROR: required --exe not defined\n\n",
					-exitval => 2); 
	}

	# make sure required directories exist
	if ( ! defined $sample_dir ) {
		pod2usage(-message => "ERROR: required --sample_dir not defined\n\n",
					-exitval => 2); 
	}
	if ( ! -d $sample_dir ) { 
		pod2usage(-message => "ERROR: --sample_dir is not a directory\n\n",
					-exitval => 2); 
	}
	
	# check split_num
	if ( ! defined $split_num ) {
		$split_num = 10000000;
		$logger->info("Setting --num to $split_num");
	}
	
	$logger->info("--exe: $exe");
	$logger->info("--sample_dir: $sample_dir");
	$logger->info("--num: $split_num");
	
	return 1;
}

sub _is_defined {
    my ($val, $default) = @_;
    
    if ( defined $val ) {
        return $val;
    }
    elsif ( defined $default ) {
        return $default;
    }
    else {
        return undef;
    }
}

sub stall {
    my ($job_name, $out_dir) = @_;
    
    $logger->info("Stalling until <$job_name> jobs finish");
    
    my $command = 'bjobs -J ' . $job_name . " > $out_dir/ACTIVE_JOBS";
    `$command`; # run the command
    
    if ( -s "$out_dir/ACTIVE_JOBS" ) {
        return 1;  # send stall singnal
    }
    
    # remove the ACTIVE_JOBS file before sending the continue signal
    `rm "$out_dir/ACTIVE_JOBS"`;
    
    return 0;  # send continue signal
}

__END__

# POD

=head1 NAME

run_rRNA_pred.pl - Runs rRNA_prediction to estimate domain level proportions


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    run_rRNA_pred.pl
        -e ~/apps/rRNA_prediction/run.sh
        -d my_metagenome_samples/
		-fi file_identifier
		[seqtk /proj/dangl_lab/apps/seqtk/seqtk]
        [-n 10000000]
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --exe | -e			Path to run_rRNA_prediction.sh 
    --sample_dir | -d	Path to root directory where samples are stored
	--file_ident | -fi	String that identifies the metagenome file to use
    --num | -n			Number of reads to put in each split sub file
	--seqtk				Path to the seqtk script
    --help | -h     	Prints USAGE statement
    --man           	Prints the man page
    --debug	        	Prints Log4perl DEBUG+ messages
    --verbose       	Prints Log4perl INFO+ messages
    --quiet	        	Suppress printing ERROR+ Log4perl messages
    --logfile       	File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --exe | -e

The exe path to the run.sh command for rRNA_Predict
    
=head2 --sample_dir | -d

Path to directory with sample data.  Each sample should have a subdirectory
inside --sample_dir.  Within each sample's directory should be a fasta file with
reads.  These fasta files will be split into many sub directories within the
sample directory.

=head2 --file_ident | -fi

This is a unique string from the metagenome files in the directories. This ensures the correct files are used for analysis.

=head2 [--seqtk]

Path to the executable seqtk. This will be used to convert fastq files into fasta files. If you have your files formatted as fasta files already this option does not have to be specified

=head2 [--num | -n]

Each sample is split into smaller files to make the jobs managable.  --num is
the number of sequences to put in each split file.
DEFAULT: 10000000
 
=head2 [--help | -h]
    
An optional parameter to print a usage statement.

=head2 [--man]

An optional parameter to print he entire man page (i.e. all documentation)

=head2 [--debug]

Prints Log4perl DEBUG+ messages.  The plus here means it prints DEBUG
level and greater messages.

=head2 [--verbose]

Prints Log4perl INFO+ messages.  The plus here means it prints INFO level
and greater messages.

=head2 [--quiet]

Suppresses print ERROR+ Log4perl messages.  The plus here means it suppresses
ERROR level and greater messages that are automatically printed.

=head2 [--logfile]

File to save Log4perl messages.  Note that messages will also be printed to
STDERR.
    

=head1 DESCRIPTION

I want to determine the rough percentage of each domain that make up my
metagenome samples (ie 60% bacterial, 30% eukaryote, etc).  Ruben (from Paul's)
lab suggested the rRNA_prediction software.  I've decided that I don't like it
for reasons that I explain on the trello card "Identification of ribosomal RNA
genes in metagenomic fragments" on my "Papers" board.  

=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Pod::Usage
Carp
Readonly
version
Log::Log4perl qw(:easy)
Log::Log4perl::CommandLine qw(:all)
File::Basename
BioUtils::FastaIO


=head1 AUTHORS

Scott Yourstone     scott.yourstone81@gmail.com
Nick Colaianni		nick.colaianni@gmail.com
    
    
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2015, Scott Yourstone
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.


=cut
