#! /usr/bin/env perl

#This program will add an extra column to the log ratio that basically just says if it has a a fxn or not.
#You need to pass the cog metadata file and the log_ratio file

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
my $cog_file;
my $log_ratio;
my $out;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'cog_file|cf=s'  =>  \$cog_file,
            'log_ratio|lf=s' =>  \$log_ratio,
            'out|o=s'        =>  \$out,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

#Main

#read in cog file
my $cog_fo = file( $cog_file );
my @cog_file = $cog_fo->slurp( split=>qr/\t/, chomp=>1 );
my $hdr = shift @cog_file;

#read lr file
my $lr_fo = file( $log_ratio );
my @lr_file = $lr_fo->slurp( split=>qr/\t/, chomp=>1 );

my $cog2fxn_href = create_cog_fxn_link( \@cog_file );

print_out_new_lr_file( \@lr_file, $cog2fxn_href );


## Subroutines ##

sub create_cog_fxn_link {
    my ( $cog_aref ) = @_;
    
    my %hash;
    foreach my $line_aref ( @$cog_aref ) {
        if ( $line_aref->[2] =~ qr/uncharacterized/i ) {
            $hash{$line_aref->[0]} = "0";
        }
        else {
            $hash{$line_aref->[0]} = "1";
        }
    }
    
    return \%hash;
}

sub print_out_new_lr_file {
    my ( $lr_aref, $href ) = @_;
    
    open my $OUT, ">", $out;
    my $head_aref = shift @$lr_aref;
    
    print $OUT join("\t", @$head_aref), "\tfxn\n";
    
    foreach my $line_aref ( @$lr_aref ) {
        if ( $href->{$line_aref->[1]} ) {
            print $OUT join("\t", @$line_aref), "\t", $href->{$line_aref->[1]}, "\n";
        }
        else {
            print $OUT join("\t", @$line_aref), "\t0\n";
        }
    }
    
    close $OUT;
    return 1;
}
