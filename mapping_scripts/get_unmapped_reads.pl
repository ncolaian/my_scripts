#! /usr/bin/env perl

# Must have samtools module uploaded

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# My Variables
my $help = 0;
my $man = 0;
my $bam_dir;
my $sample_read_dir;
my $sample_names;
my $out_dir;


# Read in variables from the command line
GetOptions( 'man'   =>  \$man,
           'help'   =>  \$help,
           'bam_dir|bd=s'       =>  \$bam_dir,
           'out_dir|o=s'        =>  \$out_dir,
           'samples|s=s'        =>  \$sample_names,
           'smpl_rd_dir|srd=s'  =>  \$sample_read_dir,
           ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
# Get the genome and sample names into arrays
$logger->info("getting sample and genome names from the files");
my $samples_list_aref = read_file_into_an_array( $sample_names );

#Go through each bam file and extract information
open my $OUT, ">", "$out_dir/unused_reads.fasta";
foreach my $smpl ( @$samples_list_aref ) {
    my $mapped_reads_href = get_mapped_reads( $smpl );
    get_unmapped_reads( $OUT, $smpl, $mapped_reads_href );
}

close( $OUT );
## Subroutines ##

sub read_file_into_an_array {
    my ($file) = @_;
    
    my $file_o = file($file);
    my @file_list = $file_o->slurp( chomp=>1 );
    
    return \@file_list;
}

sub get_mapped_reads {
    my ( $samp ) = @_;
    $logger->info( "Getting the mapped reads for $samp" );
    
    my $bam_path = "$bam_dir/$samp.bam";
    my %hash; # This hash will hold all the reads that are mapped
    
    #Open the bam file to read -- Must have samtools
    open(my $BAM, "samtools view $bam_path |") || croak "Could not use samtools to open up $bam_path";
    
    foreach my $line ( <$BAM> ) {
        chomp $line;
        my @split_ln = split /\t/, $line;
        $hash{$split_ln[3]} = 1; # Gets the ID for a mapped read
    }
    
    return \%hash;
}

sub get_unmapped_reads {
    my ( $FH, $sample_id, $mapped_href ) = @_;
    $logger->info( "Printing $sample_id\'s unmapped reads to the output fasta" );
    
    # read original read file into an array
    my $orig_read_file = read_file_into_an_array( "$sample_read_dir/$sample_id/$sample_id.raw.fastq" );
    
    #read the id's in the fastq file
    for ( my $i = 0; $i < scalar( @$orig_read_file ); $i = $i + 4 ) {
        my $id_line = $orig_read_file->[$i];
        chomp $id_line;
        next if ( $mapped_href->{$id_line} );
        
        #Print out the reads in fasta format
        my $sequence = $orig_read_file->[$i+1];
        print $FH ">$id_line\n$sequence\n";
    }
    
    return 1;
}
