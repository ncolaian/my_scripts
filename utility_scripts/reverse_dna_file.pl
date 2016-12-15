#! usr/bin/env perl

#Get reverse compliments of everything in the file, will print results in the same directory

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

# My Variables
my $help = 0;
my $man = 0;
my $file;


# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'file=s'     => \$file,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#subroutines



#Main
my $file_obj = file( $file );
my @slurped_file = $file_obj->slurp( chomp=>1 );

$file =~ s/.fasta//;
open ( my $FH, ">", "$file.rc.fasta");

foreach my $line ( @slurped_file ) {
    if ( $line =~ qr/>/ ) {
        print $FH $line, "\n";
    }
    else {
        my $rv_line = reverse( $line );
        $rv_line = uc $rv_line;
        $rv_line =~ tr/ACTG/TGAC/;
        print $FH $rv_line, "\n";
    }
}


