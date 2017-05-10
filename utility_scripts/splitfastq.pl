#! /usr/bin/perl

#This program will split a fastq file that has both of the paired end reads

use strict;
use warnings;

#Check the arguements being passed in to the program
if (scalar(@ARGV) != 3){
    print STDERR "Arguements to pass: <Fastq_file>  <1_output> <2_output>\n";
    exit(1);
}

#Read in the arguements 

my $fastq = $ARGV[0];
my $output1 = $ARGV[1];
my $output2 = $ARGV[2];

#Open Files

open(FILE,$fastq)||die("Could not open file\n");
open(OUTFILE1,'>',"$output1.fasta");
open(OUTFILE2, '>',"$output2.fasta");

my $line;
my $tf = 0;

#Read through each line and separate the forward and backwards reads
while ($line = <FILE>){
    chomp($line);
    if ($line =~ /.*\/1/){
	print OUTFILE1 "$line\n";
	$tf = 1;
    }
    elsif ($line =~ /.*\/2/){
	print OUTFILE2 "$line\n";
	$tf = 2;
    }
    elsif ($tf == 1){
	print OUTFILE1 "$line\n";
    }
    elsif ($tf == 2){
	print OUTFILE2 "$line\n";
    }
}

close FILE;
close OUTFILE1;
close OUTFILE2;


