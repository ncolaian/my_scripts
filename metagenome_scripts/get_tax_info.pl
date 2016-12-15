#! usr/bin/env perl

# This will get the taxonomy information for each genome from the file /proj/dangl_lab/data/databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt. It will produce a tabulated file with taxonomy info

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
my $gene_name_file;
my $dir_name_file;
my $dir_path;
my $out_file;


# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'gnf=s'     => \$gene_name_file,
            'dnf=s'     => \$dir_name_file,
            'dp=s'      => \$dir_path,
            'out|o=s'   => \$out_file,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#subroutines
sub set_taxonomy_hash {
    my ( $href, $file ) = @_;
    my $file_o = file($file);
    my @slurped_file = $file_o->slurp( chomp=>1, split=>qr/;/ );
    
    foreach my $tax_aref ( @slurped_file ) {
        my $genus;
        my @full_taxa;
        foreach my $tax_pos ( @$tax_aref ) {
            my @split = split qr/__/, $tax_pos;
            if ( $split[0] =~ qr/s/ ) {
                last;
            }
            my $tax = $split[-1];
            $tax =~ s/\[//;
            $tax =~ s/\]//;
            if ( $split[0] =~ qr/g/ && scalar(@split) > 1 ) {
                $genus = $tax;
            }
            push @full_taxa, $tax;
        }
        if ( defined $genus ) {
            if ( !defined $href->{$genus} ) {
                $href->{$genus} = \@full_taxa;
            }
        }
    }
    return 1;
}

sub get_taxonomy_info {
    my ( $genome, $dir, $path, $href ) = @_;
    #get the most related genus from RAST data
    my $full_path = "$path/$dir/closest.genomes";
    if ( !-e $full_path ) {
        print "$genome doesn't have a closest.genomes file";
        return;
    }
    my $close_file = file( $full_path );
    my @slurp_close = $close_file->slurp( chomp=>1, split=>qr/\t/ );
    my $first_line_aref = $slurp_close[0];
    my $name = $first_line_aref->[-1];
    my @split_name = split qr/\s/, $name;
    my $genus = $split_name[0];
    
    my $full_line;
    if ( !defined $href->{$genus} ) {
        carp "$genome\'s top genus match ($genus) is not in the taxonomy file";
    }
    else {
        my $tax_info = $href->{$genus};
        $full_line = "$genome\t" . join("\t",@$tax_info) . "\n"; 
    }
    return $full_line;
}
# main
my $taxonomy_href = {}; # Will hold all information up to genus w/ genus as key
my $tax_file = "/proj/dangl_lab/data/databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt";

set_taxonomy_hash( $taxonomy_href, $tax_file );

my $gene_nfo = file( $gene_name_file );
my $dir_nfo = file( $dir_name_file );
my @slurp_gene = $gene_nfo->slurp( chomp=> 1 );
my @slurp_dir = $dir_nfo->slurp( chomp=> 1 );

open( my $FH, ">", $out_file );
print $FH "Genome_Name\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\n";
for ( my $i = 0; $i < scalar(@slurp_gene); $i++ ) {
    my $gene = $slurp_gene[$i];
    my $dir = $slurp_dir[$i];
    
    my $info = get_taxonomy_info( $gene, $dir, $dir_path, $taxonomy_href );
    if ( defined $info ) {
        print $FH $info;
    }
}

close $FH;


