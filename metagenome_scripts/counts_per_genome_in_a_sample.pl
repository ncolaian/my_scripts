#! /usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# My Variables
my $help = 0;
my $man = 0;
my $count_dir;
my $ordered_genomes;
my $sample_ids;
my $file_name;
my $outfile;


# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'count_dir|cd=s'        => \$count_dir,
            'ordered_id_file|gi=s'  => \$ordered_genomes,
            'sample_ids|si=s'       => \$sample_ids,
            'file_name|f=s'         => \$file_name,
            'outfile|o=s'           => \$outfile,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

# Main #

$logger->info( "Getting sample ids" );
open my $SIN, "<", $sample_ids;
my @sample_order;
foreach my $sid ( <$SIN> ) {
    chomp $sid;
    push @sample_order, $sid;
}

$logger->info( "Getting ordered genome ids and getting count info" );
open my $FH, "<", $ordered_genomes;
my %counts_hash;
my @genome_order;
foreach my $id ( <$FH> ) {
    chomp $id;
    get_counts( \@sample_order, $id, \%counts_hash, $count_dir, $file_name );
    push @genome_order, $id;
}
close($FH);

$logger->info( "Printing out count information" );
open my $OUT, ">", $outfile;
print $OUT "SampleID\t", join("\t", @genome_order), "\n";

foreach my $sample ( @sample_order ) {
    print $OUT $sample;
    foreach my $gid ( @genome_order ) {
        my $count = $counts_hash{$gid}->{$sample};
        print $OUT "\t$count";
    }
    print $OUT "\n";
}
close $OUT;

# Subroutines #
sub get_counts {
    my ( $sample_order_aref, $gene_id, $count_href, $dir, $name ) = @_;
    open my $CFH, "<", "$dir/$gene_id/$name";
    my $first_line = <$CFH>;
    foreach my $line ( <$CFH> ) {
        chomp $line;
        my @split_line = split /\t/, $line;
        for (my $i=1; $i < scalar( @split_line ); $i++) {
            $count_href->{$gene_id}->{$sample_order_aref->[$i-1]} += $split_line[$i];
        }
    }
    close $CFH;
}