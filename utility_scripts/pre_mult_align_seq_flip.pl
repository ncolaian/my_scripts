#! usr/bin/env perl

#This script will make sure all the sequences to be aligned are on the same strand. This should be performed before multiple sequence alignment but after using aligners to pull out similar sequences. Will put the sequences on the same strand as the first sequence
#need to have my module Master_aln.pm

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Data::Dumper;
use Statistics::Descriptive;
use File::Temp qw/ tempfile tempdir /;
use List::Util qw/min max/;
use lib "/nas02/home/n/c/ncolaian/my_scripts/lib";
use Master_aln;


# Variables
my $help = 0;
my $man = 0;
my $alignment_file;
my $out_file;

#Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'af=s'      => \$alignment_file,
            'out_f|o=s' => \$out_file,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#Subroutines
sub reverse_seq {
    my ( $strand ) = @_;
    my $reverse = reverse( $strand );
    $reverse = uc $reverse;
    $reverse =~ tr/ACTG/TGAC/;
    return $reverse;
}

sub get_max_score {
    my ( $sequence, $orig_seq ) = @_;
    my $reverse_seq = reverse_seq( $sequence );
    my @split_orig_seq = split //,$orig_seq;
    my @split_sequence = split //,$sequence;
    my @split_reverse_seq = split //,$reverse_seq;
    my $orig_strand_obj = Master_aln->new( ( 'dnaseqi_ref' => \@split_orig_seq,
                                             'dnaseqj_ref' => \@split_sequence,
                                             'match_score' => 1,
                                             'mismatch_score' => -1,
                                             'gap_score'   => -2,) );
    my $reverse_strand_obj = Master_aln->new( ( 'dnaseqi_ref' => \@split_orig_seq,
                                                'dnaseqj_ref' => \@split_reverse_seq,
                                                'match_score' => 1,
                                                'mismatch_score' => -1,
                                                'gap_score'   => -2,) );
    my ( $throw_away, $pointless, $orig_max ) =
    $orig_strand_obj->get_max_coordinates();
    
    my ( $throw_away2, $pointless2, $reverse_max ) =
    $reverse_strand_obj->get_max_coordinates();
    
    if ( $orig_max >= $reverse_max ) {
        return $sequence;
    }
    else{
        return $reverse_seq;
    }
}

sub print_alignment {
    my ( $header, $seq, $FH ) = @_;
    print $FH "$header\n$seq\n";
}

#main
my $align_fo = file( $alignment_file );
my @align_slurp = $align_fo->slurp( chomp=>1 );
open( my $OFH, ">", $out_file );

print_alignment( $align_slurp[0], $align_slurp[1], $OFH );

for ( my $i = 3; $i < scalar( @align_slurp ); $i = $i + 2 ) {
    my $max_seq = get_max_score( $align_slurp[$i], $align_slurp[1] );
    print_alignment( $align_slurp[$i-1], $max_seq, $OFH );
}




