#! /usr/bin/env perl

#Combine poolcount files together into one file

use strict;
use warnings;

use Carp;
use Path::Class;
use Data::Dumper;

my ( $file1, $file2, $output_file) = @ARGV if scalar(@ARGV) == 3 || die "Did not pass in 2 poolcount files and an output file\n";

open my $FH1, "<", $file1;
open my $FH2, "<", $file2;
open my $OUT, ">", $output_file;

my $h1 = <$FH1>;
my $h2 = <$FH2>;
chomp $h1;
chomp $h2;

my @header = split /\t/, $h1;
my @add_header = split /\t/, $h2;

#push the header without checking the first column
for ( my $i = 5; $i < scalar(@add_header); $i++) {
    push @header, $add_header[$i];
}

print $OUT join( "\t", @header ), "\n";

while ( my $line1 = <$FH1> ) {
    my $line2 = <$FH2>;
    chomp $line1;
    chomp $line2;
    
    my @line1_ar = split /\t/, $line1;
    my @line2_ar = split /\t/, $line2;
    
    if ( $line1_ar[0] ne $line2_ar[0] ) {
        die "The order of the pool files do not match";
    }
    
    for ( my $c = 5; $c < scalar( @line2_ar ); $c++ ) {
        push @line1_ar, $line2_ar[$c];
    }
    
    #print out line
    print $OUT join( "\t", @line1_ar ), "\n";
}

close($FH1);
close($FH2);
close( $OUT );

