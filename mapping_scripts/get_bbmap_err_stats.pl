#! /usr/bin/env perl

# Script creates a table of 

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
my $db_file;
my $as_file;
my $out_file;

# Read in variables from the command line
GetOptions( 'man'   =>  \$man,
           'help'   =>  \$help,
           'db_file|db=s'   =>  \$db_file,
           'as_file|as=s'   =>  \$as_file,
           'out_file|o=s'   =>  \$out_file,
           ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
open my $OPEN, "<", $as_file || die ( "Cannot read $as_file\n" );
open my $OUT, ">", $out_file || die ( "Cannot open $out_file for writing\n" );

# Print header of table
print $OUT "SampleID\tPercent\tTotal\tMap_type\n";

# get assembly info and total counts
my %total;

my $position = 1;
my $genome;
foreach my $line ( <$OPEN> ) {
  chomp $line;
  if ($position == 1) {
    print $OUT "$line\t";
    $genome = $line;
    $position++;
  }
  
  elsif ( $position == 2 )  {
    my @split = split /\(/, $line;
    #percent
    my @sec_split = split /:/, $split[1];
    my $perc = $sec_split[0];
    chop $perc;
    chop $perc;
    print $OUT "$perc\t";

    #total reads
    my @tot_split = split /\+/, $line;
    chop $tot_split[0];
    print $OUT $tot_split[0], "\tAS\n";
    $position++;
  }
  
  else{
    my @split = split /\+/, $line;
    chomp $split[0];
    $total{$genome} = $split[0];
    $position = 1;
  }
}
close( $OPEN );

open my $FH, "<", $db_file || ( "Cannot read $db_file\n" );

#get db
my $sample;
my $pct_rds;
my $num_reads;
#my $pct_col = 0;
foreach my $line ( <$FH> ) {
    #handles sample name
    if ( $line =~ m/^samtools/ ) {
        print "hi";
        my @split_line = split /\./, $line;
        my @spliter = split /\s/, $split_line[0];
        $sample = $spliter[scalar(@spliter)-1];
    }
    
    #finds column with pct reads
    #my $col =0;
    #if ( $line =~ qr/^Read\s1\sdata/ ) {
    #    my @split_line = split /\t/, $line;
    #    foreach my $split ( @split_line ) {
    #        if ( $split =~ qr/^pct\sreads$/ ) {
    #            $pct_col = $col;
    #            last;
    #        }
    #        $col++;
    #    }
    #}
    
    #finds number of mapped reads
    if ( $line =~ qr/mapped\s\(100/ ) {
        my @split_line = split /\+/, $line;
        chop( $split_line[0] );
        my $full_reads = $total{$genome};
        $num_reads = $split_line[0];
        $pct_rds = ($num_reads/$full_reads)*100;
        
        #prints out info at last line
        print $OUT "$sample\t$pct_rds\t$num_reads\tDB\n";
    }
}

close ( $OUT );
close ( $FH );
