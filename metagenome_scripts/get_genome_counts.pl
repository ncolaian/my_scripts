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
my $genome_names;
my $file_name;
my $outfile;


# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'count_dir|cd=s'    => \$count_dir,
            'genome_names|gn=s' => \$genome_names,
            'file_name|f=s'       => \$file_name,
            'outfile|o=s'       => \$outfile,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

# Main #
my %genome_counts;

$logger->info( "counting up the counts in all the files" );
open my $GNIN, "<", $genome_names;
foreach my $name ( <$GNIN> ) {
    chomp $name;
    $genome_counts{$name} = 0;
}
close($GNIN);

get_counts_from_count_file( \%genome_counts, $file_name );

print_count_info( \%genome_counts );

## Subroutines ##

sub get_counts_from_count_file {
    my ( $count_href, $file ) = @_;
    foreach my $genome ( keys %$count_href ) {
        $logger->info("Getting counts from $genome");
        open (CIN, "<", "$count_dir/$genome/$file") || die "File doesn't exist at: $count_dir/$genome/$file_name";
        
        foreach my $line ( <CIN> ) {
            next if ( $line =~ qr/^__/ || $line =~ qr/^name/i);
            chomp $line;
            my @split_line = split /\t/, $line;
            for ( my $i = 1; $i < scalar( @split_line ); $i++ ) {
                my $count = $count_href->{$genome};
                $count_href->{$genome} = $count + int($split_line[$i]);
            }
        }
        close(CIN);
    }
    return 1;
}

sub print_count_info {
    my ( $count_href ) = @_;
    $logger->info( "Printing out total count file" );
    open my $OUT, ">", $outfile;
    print $OUT "Genome Id\tTotal Counts\n";
    
    foreach my $genome ( keys %$count_href ) {
        my $count = $count_href->{$genome};
        print $OUT "$genome\t$count\n";
    }
    close($OUT);
    return 1;
}



