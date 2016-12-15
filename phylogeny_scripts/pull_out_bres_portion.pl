#! usr/bin/env perl

#Tab file must follow this pattern: QueryId, SbjID, Qstart, Qend, Score, Total Length, # of Identical matches

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
my $blast_tab_file;
my $assembly_file;
my $out;
my $name;


# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'btf=s'     => \$blast_tab_file,
            'af=s'      => \$assembly_file,
            'out|o=s'   => \$out,
            'name=s'    => \$name,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#subroutines

sub print_out_blast_seq {
    my ( $sequence, $id, $FH, $id_href, $assembly_file, $name ) = @_;
    my $start_end_aref = $id_href->{$id};
    my $start = $start_end_aref->[0];
    my $end = $start_end_aref->[1];
    my @full_sequence = split //, $sequence;
    
    print $FH ">$name start=$start;end=$end;\n";
    
    for ( my $i = $start-1; $i < $end; $i++) {
        print $FH $full_sequence[$i];
    }
    print $FH "\n";
}


# Main

my %query_ids;

#Go through blast results and create a hash of ids found
my $tab_file_obj = file( $blast_tab_file );
my @divided_lines = $tab_file_obj->slurp( chomp=>1, split=>qr/\t/ );

foreach my $line_aref ( @divided_lines ) {
    if ( $line_aref->[2] < $line_aref->[3] ) {
        $query_ids{$line_aref->[0]} = [$line_aref->[2], $line_aref->[3]];
    }
    else {
        $query_ids{$line_aref->[0]} = [$line_aref->[3], $line_aref->[2]];
    }
}
#Go through assembly file and pull out the full sequences
open ( my $BSFILE, ">", $out );
my $ass_fo = file( $assembly_file );
my @ass_lines = $ass_fo->slurp( chomp=>1 );
my $t_or_f = 0;
my @sequence;
my $id;

foreach my $line ( @ass_lines ) {
    if ( $line =~ qr/>/ ) {
        my @mm_array = split qr/\s/, $line;
        if ( $t_or_f == 1 ) {
            $t_or_f = 0;
            my $full_seq = join "", @sequence;
            print_out_blast_seq( $full_seq, $id, $BSFILE, \%query_ids, $assembly_file, $name);
            undef @sequence;
        }
        
        $id = $mm_array[0];
        $id =~ s/>//;
        if ( exists $query_ids{$id} ) {
            $t_or_f = 1;
        }
    }
    
    elsif ( $t_or_f == 1 ) {
        push @sequence, $line;
    }
}

if ( $t_or_f == 1 ) {
    my $full_seq = join "", @sequence;
    print_out_blast_seq( $full_seq, $id, $BSFILE, \%query_ids, $assembly_file, $name );
}


close( $BSFILE );

