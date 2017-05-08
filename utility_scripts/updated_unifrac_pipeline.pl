#! /usr/bin/env perl

#This is the unifrac pipeline. Originally created by Scott Yourstone
#Requires qiime/1.5.0
#Made adjustments to the pipeline for full analysis in one script
#Now just need to pass in a fasta file

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Path::Class;
use Data::Dumper;
use MyX::Generic;
use Bio::SeqIO;
use File::Temp qw(tempfile);
use File::Basename;
use Scalar::Util qw(looks_like_number);
use File::Copy;

# my subroutines
sub stall;
sub make_otus;
sub make_otu_table;
sub filter_contam;
sub filter_otu_reads;
sub filter_otu_samples;
sub check_OTU_headers;
sub make_filt_rep_seqs;
sub convert_to_biom;
sub build_msa;
sub build_phylogeny;
sub run_unifrac;
sub calculate_pcoa;
sub reformat_pcoa_files;
sub run_alpha_diversity_command;


#my Variables
my $help;
my $man;
my $fasta_file;
my $otu_dir;
my $unifrac_dir;
my $min_read_filt = 0;
my $min_sample_filt = 0;
my $make_otu_prg_path;

# Read in the variables
GetOptions ( 'man' => \$man,
             'help' => \$help,
             'uc_fa|uf=s'   => \$fasta_file,
             'otu_dir|od=s' => \$otu_dir,
             'uni_dir|ud=s' => \$unifrac_dir,
             'min_reads=i'    => \$min_read_filt,
             'min_samples=i'  => \$min_sample_filt,
             'make_otu_path|mop=s' => \$make_otu_prg_path,
             ) || die("There was an error in the command line arguements");

# Pod Usage for the Manual and Help pages.
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage( -verbose => 3 ) }

# Setup logging environment
my $logger = get_logger();

#Check and handle the filter params
if ( $min_read_filt != 0 || $min_sample_filt != 0 ) {
    $logger->info("Checking the filter params");
    if ( $min_read_filt == 0 || $min_sample_filt == 0 ) {
        croak( "You must pass in both the minimum samples and minimum read parameters if you want to filter the tables before analysis" );
    }
    if ( ! looks_like_number( $min_read_filt ) && ! looks_like_number( $min_sample_filt ) ) {
        croak( "The filter paramaters passed in must be in numeric integer form. i.e 1" );
    }
}

## SUBROUTINES ##
sub stall {
    my ($job_name, $out_dir) = @_;
    
    $logger->info("Stalling until <$job_name> jobs finish");
    
    my $command = 'bjobs -J ' . $job_name . " > $out_dir/ACTIVE_JOBS";
    `$command`; # run the command
    
    if ( -s "$out_dir/ACTIVE_JOBS" ) {
        return 1;  # send stall signal
    }
    
    # remove the ACTIVE_JOBS file before sending the continue signal
    `rm "$out_dir/ACTIVE_JOBS"`;
    
    return 0;  # send continue signal
}

sub make_otus {
    $logger->info( "Making OTU's with fasta file" );
    my $cmd = "bsub -J 1.2 -M 35 -o $otu_dir/lsf.out " . "perl $make_otu_prg_path/make_otus.pl -i $fasta_file -denovo_chimera -ref_chimera -minsize 2 -M uparse -o $otu_dir";
    $logger->debug($cmd);
    system($cmd);
}

sub make_otu_table {
    $logger->info( "Making an OTU table with newly created OTU's" );
    my $cmd = "perl $make_otu_prg_path/make_otu_table.pl -i $otu_dir/otu_table.uc -s $otu_dir/rep_set.fa -o $otu_dir/otu_table";
    $logger->debug($cmd);
    system($cmd);
}

sub filter_contam {
    $logger->info ( "Filter contaminants" );
    my $cmd = "perl /proj/dangl_lab/bin/filter_contaminants.pl -blast_db /nas02/home/y/o/yourston/dangl_lab/data/contaminants_db/contaminants.fasta -query_file $otu_dir/rep_set.fa -output_dir $otu_dir -otu_table $otu_dir/otu_table/otu_table.txt";
    $logger->debug($cmd);
    system($cmd);
}

sub filter_out_min {
    $logger->info( "Filtering out the OTU's that don't meet minimum requirements" );
    open my $FH, "<", "$otu_dir/out_non_contaminants_otu_table.txt";
    open my $OUT, ">", "$otu_dir/otu_table_filt.txt";
    
    my $skip_line = <$FH>;
    my $first_line = <$FH>;
    print $OUT $skip_line;
    print $OUT $first_line;
    
    my $amn_filtered = 0;
    
    foreach my $line ( <$FH> ) {
        chomp $line;
        #keep track of the reads that are below threshold
        my $pass = 0;
        $line =~ s/\.0//g;
        my @line_a = split /\t/, $line;
        #go through and count reads that are above min
        foreach my $count ( @line_a ) {
            if ( $count =~ /^[0-9]*$/m ) {
                if ( $count >= $min_read_filt ) {
                    $pass++;
                }
            }
        }
        
        #Only print the lines that pass the filter
        if ( $pass >= $min_sample_filt ) {
            print $OUT $line, "\n";
        }
        else {
            $amn_filtered++;
        }
    }
    print "Amount of OTUS filtered: $amn_filtered\n";
    close( $FH );
    close( $OUT );
}

sub check_OTU_headers {
    my ($file) = shift @_; 

    # open the OTU table file
    open my $IN, "<", $file or die $!; 
    
    my $skip_line = <$IN>;
    my $first_line = <$IN>;

    # if the header line already starts with OTU_ID return
    if ( $first_line =~ m/^OTU_ID/ ) { 
        return 1;
    }   
    
    # if I get to here I need to add OTU_ID\t to the headers
    # first open a tmp file
    my ($fh, $filename) = tempfile();
    my @f_line = split /\t/,$first_line;
    shift @f_line;
    print $fh "OTU_ID\t" . join "\t",@f_line;
    foreach my $line ( <$IN> ) {
        $line =~ s/\.0//g;
        print $fh $line;
    }

    close($fh);
    close($IN);

    # replace the updated file -- update the symbolic link if neccessary
	my $file_to_update = $file;
	if ( -l $file ) {
		my $link_dir = dirname($file);
		my $link = readlink($file);
		$file_to_update = "$link_dir/$link";
	}
    `mv $filename $file_to_update`;
}

sub make_filt_rep_seqs {
	$logger->info("make filtered rep seqs fasta file");
	
	my $awk = "awk \'{print \$1}\' " .
				"$otu_dir/otu_table_filt.txt > " .
				"$otu_dir/otu_table_filt_names.txt";
	$logger->debug($awk);
	`$awk`;
	
	my $get_seqs = "perl /nas02/home/y/o/yourston/scripts/SeqTools/fasta_get_seqs_by_id.pl " .
					"$otu_dir/rep_set.fa " .
					"$otu_dir/otu_table_filt_names.txt " .
					"$otu_dir/rep_set_filt.fa";
	$logger->debug($get_seqs);
	`$get_seqs`;
}

sub convert_to_biom {
	$logger->info("Convert to otu table back to biom format");
	
	if ( -s "$otu_dir/otu_table_filt.biom" ) {
		`rm $otu_dir/otu_table_filt.biom`;
	}
	
	my $cmd = "biom convert " .
				"-i $otu_dir/otu_table_filt.txt " .
				"-o $otu_dir/otu_table_filt.biom " .
				"--table-type=\"otu table\"";
	$logger->debug($cmd);
	`$cmd`;
}

sub build_msa {
	$logger->info("Build MSA");
	
	my $cmd = "align_seqs.py " .
				"-i $otu_dir/rep_set_filt.fa " .
				"-o $unifrac_dir";
	$logger->debug($cmd);
	`$cmd`;
}

sub build_phylogeny {
	$logger->info("Build phylogeny");

	my $cmd = "make_phylogeny.py " .
				"-i $unifrac_dir/rep_set_filt_aligned.fasta";
	$logger->debug($cmd);
	`$cmd`;
    system( "mv taxonomy/* $unifrac_dir");
    system( "rm -r taxonomy");
}

sub run_unifrac {
	$logger->info("run unifrac");
	
	my $cmd = "beta_diversity.py " .
				"-i $otu_dir/otu_table_filt.biom " .
				"-o $unifrac_dir/ " .
				"-t $unifrac_dir/rep_set_filt_aligned.tre";
	$logger->debug($cmd);
	`$cmd`;
}

sub calculate_pcoa {
	# weighted
	$logger->info("Cacluate weighted pcoa");
	
	my $cmd = "principal_coordinates.py " .
				"-i $unifrac_dir/weighted_unifrac_otu_table_filt.txt " .
				"-o $unifrac_dir/tmp_weighted_pcoa.txt";
	$logger->debug($cmd);
	`$cmd`;

	# unweighted
	$logger->info("Cacluate unweighted pcoa");

	$cmd = "principal_coordinates.py " .
			"-i $unifrac_dir/unweighted_unifrac_otu_table_filt.txt " .
			"-o $unifrac_dir/tmp_unweighted_pcoa.txt";
	$logger->debug($cmd);
	`$cmd`;
}

sub reformat_pcoa_files {
	my ($tmp_file, $pcoa_file, $eigen_file) = @_;

	# some edits to the pcoa file
	`sed -i 's/ vector number//g' $tmp_file`;
	
	# seperate the last 2 lines of the pcoa file into a seperate file
	open my $IN, "<", "$tmp_file" or
		$logger->logdie("Cannot open $tmp_file");
	
	my @lines = <$IN>;
	
	open my $EGN, ">", "$eigen_file" or
		$logger->logdie("Cannot open file: $eigen_file");
	
	# print the eigen values and percent variation explained
	print $EGN (pop @lines), "\n";
	print $EGN (pop @lines);
	
	# remove the white spaces
	pop @lines;
	pop @lines;
	
	# print the pcoa data
	open my $PCOA, ">", $pcoa_file or
		$logger->logdie("Cannot open file: $pcoa_file");
	print $PCOA @lines;
	
	# remove the tmp file
	`rm $tmp_file`;
}

sub run_alpha_diversity_command {
	my $cmd;
	
	$logger->info("Run alpha diversity command");
	$cmd = "alpha_diversity.py " .
		"-i $otu_dir/otu_table_filt.biom " .
		"-o $unifrac_dir/alpha_diversity.txt " .
		"-m PD_whole_tree,chao1,shannon,simpson " .
		"-t $unifrac_dir/rep_set_filt_aligned.tre";
	
	$logger->debug($cmd);
	`$cmd`;
}

## MAIN ##

# Make sure the file passed in is a fasta file
$logger->info("checking if file is a fasta file");
my $seqio = Bio::SeqIO->new(-file => $fasta_file , -format => "fasta");

# Use the fasta file to make the otus
make_otus();
while( stall( "1.2", $otu_dir) ) {
    $logger->info( "Stalling" );
    sleep 720;
}
# Make the OTU table
make_otu_table();

#filter contaminants out of the otu table
filter_contam();

#One more manual filter step
if ( $min_read_filt == 0 && $min_sample_filt == 0 ) {
    my $cmd = "cp $otu_dir/out_non_contaminants_otu_table.txt $otu_dir/otu_table_filt.txt";
    system( $cmd );
}
else{
    $logger->info( "Filtering by the minimum requirements" );
    #produces the filtered otu table by the filter requirements passed
    filter_out_min();    
}

# check the headers to make sure there is a col named OTU_ID at the beginning
$logger->info( "Checking Header" );
check_OTU_headers("$otu_dir/otu_table_filt.txt");

#get only the OTUs that are in the filtered OTU table
make_filt_rep_seqs();

# convert the filtered OTU table back to biom format
convert_to_biom();

# build an MSA of the OTU sequences
build_msa();

# build phylogeny
build_phylogeny();

# run unifrac
run_unifrac();

# calculate pcoa
calculate_pcoa();

# reformat the pcoa file
reformat_pcoa_files(
	"$unifrac_dir/tmp_weighted_pcoa.txt",
	"$unifrac_dir/weighted_unifrac_otu_table_filt_pcoa.txt",
	"$unifrac_dir/weighted_unifrac_otu_table_filt_eigen_and_perc_var.txt"
);
reformat_pcoa_files(
	"$unifrac_dir/tmp_unweighted_pcoa.txt",
	"$unifrac_dir/unweighted_unifrac_otu_table_filt_pcoa.txt",
	"$unifrac_dir/unweighted_unifrac_otu_table_filt_eigen_and_perc_var.txt"
);


# Run the alpha diversity
run_alpha_diversity_command();

__END__

=head1 NAME

unifrac_pipeline.pl - runs unifrac

=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    unifrac_pipeline.pl
        --otu_dir|od otu_dir/
        --uni_dir|ud unifrac/
        --uc_fa|uf reads.fasta
        --make_otu_path|mop otu_programs/
        [--min_reads 3] 
        [--min samples 10] 
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]


=head1 ARGUMENTS
    
=head2 --otu_dir | -od

Path to dir where output from anything to do with otu creation will be placed. 
    
=head2 --uni_dir | -ud 

Path to where unifrac output will be stored

=head2 --uc_fa | -uf

Path to a fasta file with the data needed to create OTU's and the data that unifrac is desired to be run on.

=head2 [--min_reads]

The minimum number of counts an otu must have to be counted as significant. If min_reads is passed without a number it will not filter anything but otus with 0 counts. 

=head2 [--min_samples]

The minimum number of significant counts that an otu must have to be retained for downstream analysis. If min_reads is passed without a number 0 will be used. Meaning if there is one significant count the otu will be retained.

** If min_reads and min_samples are not passed no filtering will occur **

=head2 --make_otu_path | -mop

Path to dir where the programs for OTU creation are located. These programs are make_otus.pl and make_otu_tble.pl.
 
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

This runs the Unifrac pipeline starting with a fasta file.  Here are the steps:

    - create otus
    - create an otu table and rep_set fasta file
	- create filtered rep_set fasta file
    - filter out OTU's based on minimum count and sample information
	- convert the OTU table from txt format to biom
	- make the MSA
	- make the phylogeny
	- run unifrac
	- calculate the principle coordinates
	- do some file clean up and reformating

=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);
Path::Class;
Data::Dumper;
MyX::Generic;
Bio::SeqIO;
File::Temp qw(tempfile);
File::Basename;
Scalar::Util qw(looks_like_number);
File::Copy;


=head1 AUTHOR(S)

Scott Yourstone     scott.yourstone81@gmail.com
    &
Nicholas Colaianni  nick.colaianni@gmail.com
    
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