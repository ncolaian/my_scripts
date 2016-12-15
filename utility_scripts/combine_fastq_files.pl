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
my $file_names;
my $out_dir;

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'txt_file|txt=s' => \$file_names,
            'out|o=s'   => \$out_dir,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#Subs

#Main
my $file_obj = file($file_names);
my @lines = $file_obj->slurp( chomp=>1, split=>qr/=/);
my %names;
my $name;

foreach my $line_aref ( @lines ) {
    if ( !defined $line_aref->[0] ) {
        next;
    }
    if ( $line_aref->[0] eq "name" ) {
        $name = $line_aref->[1];
        $names{$name} = {};
    }
    elsif ( $line_aref->[0] =~ qr/f*_*/ ) {
        my $hash_ref = $names{$name};
        $hash_ref->{$line_aref->[0]} = $line_aref->[1];
    }
}

foreach my $nms_frm_hash ( keys %names ) {
    my $cmd1_part = "bsub -o $out_dir/lsf1.out 'cat ";
    my $cmd2_part = "bsub -o $out_dir/lsf2.out 'cat ";
    foreach my $file_type ( keys %{ $names{$nms_frm_hash} } ) {
        my $hash_ref = $names{$name};
        my $file_name = $hash_ref->{$file_type};
        if ($file_type  =~ qr/f*_1/) {          
            $cmd1_part .= "$file_name ";
        }
        if ($file_type =~ qr/f*_2/) {
            $cmd2_part .= "$file_name ";
        }
    }
    $cmd1_part .= "> $out_dir/$nms_frm_hash.pt1.fq';";
    $cmd2_part .= "> $out_dir/$nms_frm_hash.pt2.fq';";
    
    system($cmd1_part);
    system($cmd2_part);
}

