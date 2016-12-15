#! usr/bin/evn perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Data::Dumper;

# My Variables
my $help = 0;
my $man = 0;
my $fasta;
my $out;

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'fasta|f=s' => \$fasta,
            'out|o=s'   => \$out,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#Subs

#Main
my $file_obj = file($fasta);
my @lines_of_fasta = $file_obj->slurp( chomp=>1 );
my %hash_of_files;
my $previous_name;

for ( my $i = 0; $i < scalar( @lines_of_fasta ); $i++ ) {
    #handle first
    if ( $i == 0 ) {
        my $name = $lines_of_fasta[$i];
        $name =~ s/>//;
        $name .= ".fasta";
        $hash_of_files{$name} = [$i];
        $previous_name = $name;
    }
    #handle end of file
    elsif ( $i == (scalar( @lines_of_fasta ) - 1) ) {
         push( @{$hash_of_files{$previous_name}}, $i );
    }
    #handle the normal stuff
    elsif ( $lines_of_fasta[$i] =~ qr/>/ ) {
        my $name = $lines_of_fasta[$i];
        $name =~ s/>//;
        $name .= ".fasta";
        $hash_of_files{$name} = [$i];
        push( @{$hash_of_files{$previous_name}}, ($i -1) );
        $previous_name = $name;
    }
}


foreach my $file_names (keys %hash_of_files) {
    my $start_and_stop_aref = $hash_of_files{$file_names};
    my $start = $start_and_stop_aref->[0] + 1;
    my $end = $start_and_stop_aref->[1] + 1;
    my $cmd = "sed -n $start,$end" . "p $fasta > $out/$file_names";
    system($cmd);
}


