#! /usr/bin/env perl

use strict;
use warnings;

my ( $file1, $file2, $output_file) = @ARGV if scalar(@ARGV) == 3 || die "Did not pass in 2 poolcount files and an output file\n";

open my $FH1, "<", $file1;
open my $FH2, "<", $file2;
open my $OUT, ">", $output_file;

my $header1 = <$FH1>;
my $header2 = <$FH2>;
chomp $header1;
chomp $header2;

my @full_header = split /\t/, $header1;
my @adding_header = split /\t/, $header2;

if ( scalar( @full_header ) != scalar ( @adding_header ) ) {
    die "The files do not contain the same amount of columns so the merging can not work";
}

print $OUT $header1, "\n";

while ( my $line1 = <$FH1> ) {
    my $line2 = <$FH2>;
    chomp $line1;
    chomp $line2;
    
    my @line1_ar = split /\t/, $line1;
    my @line2_ar = split /\t/, $line2;
    
    if ( $line1_ar[0] ne $line2_ar[0] ) {
        die "The order of the pool files do not match";
    }
    
    my %adding_hash;
    #hold onto the value of the adding file
    for ( my $i = 5; $i < scalar( @line2_ar ); $i++ ) {
        $adding_hash{$adding_header[$i]} = $line2_ar[$i];
    }
    
    #put the largest value in the first array
    for ( my $i = 5; $i < scalar( @line1_ar ); $i++ ) {
        if ( $line1_ar[$i] < $adding_hash{$full_header[$i]} ) {
            $line1_ar[$i] = $adding_hash{$full_header[$i]};
        }
    }
    
    #print out line
    print $OUT join( "\t", @line1_ar ), "\n";
}

close($FH1);
close($FH2);
close( $OUT );
