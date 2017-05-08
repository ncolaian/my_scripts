#! /usr/bin/env perl

#This program will add 

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;

#My Variables
my $help = 0;
my $man = 0;
my $gene_file;
my $count_file;
my $out_file;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'gene_file|gf=s'    => \$gene_file,
            'count_file|cf=s'   => \$count_file,
            'out_file|o=s'      => \$out_file,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##

my $gfo = file( $gene_file );
my @gene_file_array = $gfo->slurp( split=>qr/\t/,chomp=>1 );
my $discard_hder_ln = shift @gene_file_array;

#Create a hash to search for the genes
my $gene_href = create_gene_hash( \@gene_file_array );

open my $CF, "<", $count_file;
my $header = <$CF>;
chomp $header;

open my $OUT, ">", $out_file;

print $OUT "$header\tgeneid\tgenename\n";

my $no_gene = 0;

foreach my $line ( <$CF> ) {
    chomp $line;
    my ($gene_aref, $no_gene_count) = get_gene_id_and_name($line, $gene_href);
    $no_gene += $no_gene_count;
   # next if ( $no_gene_count == 1 );
    my $gene_id = $gene_aref->[0];
    my $gene_name = $gene_aref->[1];
    
    print_to_out($OUT, $line, $gene_id, $gene_name);
}
print "Number of codes not in a gene: $no_gene\n";
close $OUT;
close $CF;


## Subs ##

sub create_gene_hash {
    my ( $g_aref ) = @_;
    
    my %href;
    
    foreach my $aref ( @$g_aref ) {
        my $scaffold = $aref->[3];
        my $start = $aref->[4];
        my $end = $aref->[5];
        my $com = "$start-$end";
        
        if ( $href{ $scaffold } ) {
            my $cur_href = $href{ $scaffold };
            $cur_href->{$com} = [$aref->[0], $aref->[8]];
        }
        else {
            my %scaf_href;
            $scaf_href{$com} = [$aref->[0], $aref->[9]];
            $href{$scaffold} = \%scaf_href;
        }
    }
    
    return \%href;
}

sub get_gene_id_and_name {
    my ( $line, $href ) = @_;
    my @line_array = split /\t/, $line;
    my $scaff = $line_array[2];
    
    if ( !defined $line_array[4] ) {
        return ["pastEnd", "No-Name"], 1;
    }
    
    my $scaf_href = $href->{$scaff};
    
    foreach my $key ( %$scaf_href ) {
        my @mm_array = split /-/, $key;
        my $start = $mm_array[0];
        my $end = $mm_array[1];
        if ( ref($start) eq "ARRAY" ) {
            $end = $start->[1];
            $start = $start->[0];
            print $scaff, "\n";
        }
        
        if ( $line_array[4] >= $start && $line_array[4] <= $end ) {
            my @id_name_array = @{$scaf_href->{$key}};
            return \@id_name_array, 0;
        }
    }
    
    return ["No-Gene","No-Name"], 1;
}

sub print_to_out {
    my ( $fh, $line, $id, $name ) = @_;
    my $new_line = "$line\t$id\t$name\n";
    print $fh $new_line;
    return 1;
}

