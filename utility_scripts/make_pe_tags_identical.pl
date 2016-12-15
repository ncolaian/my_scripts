#! usr/bin/evn perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Data::Dumper;
use Readonly;

# My Variables
my $help = 0;
my $man = 0;
my $name;
my $dir;

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'name=s' => \$name,
            'dir=s'  => \$dir,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# MAIN

Readonly::Hash my %NAMES => map { $_ => 1 } qw(
    split_pe_aa.1
    split_pe_aa.2
    split_pe_ab.1
    split_pe_ab.2
);

foreach my $f_name ( keys %NAMES ) {
    my $full_f_path = "$dir/$name/$f_name";
    open my $in_fh, "<", $full_f_path;
    open my $out_fh, ">", "$full_f_path.new";
    
    while (my $line = <$in_fh> ) {
        if ( $line =~ qr/@/ ) {
            my @part = split qr/\s/, $line;
            chop $part[0];
            chop $part[0];
            print $out_fh $part[0], "\n";
        }
        else {
            print $out_fh $line;
        }
    }
    close( $in_fh );
    close( $out_fh );
}


