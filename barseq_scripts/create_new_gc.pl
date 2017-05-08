#! /usr/bin/env perl

# This program will take the genes.GC file that has gene information and return 2 files. The two gene files that are created replace the locusId gene names with 1) Cog info and 2) Cog Functional annotations

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
my $gene_file;
my $annote_file;
my $cog_file;
my $base_out;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'gene_file|gf=s'    => \$gene_file,
            'annote_file|af=s'  => \$annote_file,
            'cog_file|cf=s'     =>  \$cog_file,
            'base_out|o=s'       => \$base_out,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##

my $gfo = file( $gene_file );
my @gene_file_array = $gfo->slurp( split=>qr/\t/,chomp=>1 );

my $annote_file_o = file( $annote_file );
my @annote_file_array = $annote_file_o->slurp( split=>qr/\t/, chomp=>1 );

my $cog_fo = file( $cog_file );
my @cog_file = $cog_fo->slurp( split=>qr/\t/, chomp=>1 );

# Need to create a hash that links genes to cogs
my $genes_to_cogs_href = create_gene_2_cog_link(\@annote_file_array);

# Need to create a hash that links cogs to annotation
my $ann_to_cogs_href = create_anns_2_cogs_link(\@cog_file);

# Create the gene file
create_the_two_files($genes_to_cogs_href, $ann_to_cogs_href, \@gene_file_array);


## Subroutines ##

sub create_gene_2_cog_link {
    my ( $annote_aref ) = @_;
    
    my $first_line_aref = shift @$annote_aref;
    #Find the locus tag and cog column
    my $locus_col = 0;
    my $cog_col = 0;
    my $count = 0;
    
    foreach my $col ( @$first_line_aref ) {
        if ( $col =~ qr/locus_tag/i) {
            $locus_col = $count;
        }
        elsif ( $col =~ qr/cog/i ) {
            $cog_col = $count;
        }
        $count++;
    }
    
    #Go through and make hash link
    my %hash;
    foreach my $line_aref ( @$annote_aref ) {
        $hash{$line_aref->[$locus_col]} = $line_aref->[$cog_col];
    }
    
    return \%hash;
}

sub create_anns_2_cogs_link {
    my ( $cog_file ) = @_;
    
    my $first_line_aref = shift @$cog_file;
    #find the cog and func columns
    my $func_col = 0;
    my $cog_col = 0;
    my $count = 0;
    
    foreach my $col ( @$first_line_aref ) {
        if ( $col =~ qr/func/i) {
            $func_col = $count;
        }
        elsif ( $col =~ qr/cog/i ) {
            $cog_col = $count;
        }
        $count++;
    }
    
    #Go through and make hash link
    my %hash;
    foreach my $line_aref ( @$cog_file ) {
        $hash{$line_aref->[$cog_col]} = $line_aref->[$func_col];
    }
    
    return \%hash;
}

sub create_the_two_files {
    my ( $g2c_href, $c2f_href, $gfa_aref ) = @_;
    
    my $first_line_aref = shift @$gfa_aref;
    #this assumes the locus id is the first column
    open my $COGOUT, ">", "$base_out\_cog_genes.GC";
    open my $FNOUT, ">", "$base_out\_func_genes.GC";
    
    print $COGOUT join("\t", @$first_line_aref), "\n";
    print $FNOUT join("\t", @$first_line_aref), "\n";
    
    foreach my $line_aref ( @$gfa_aref ) {
        #cog file
        if ( $g2c_href->{$line_aref->[0]} && $g2c_href->{$line_aref->[0]} ne "NA" ) {
            my $cog = $g2c_href->{$line_aref->[0]};
            $line_aref->[0] = $cog;
            $line_aref->[1] = $cog;
            print $COGOUT join("\t", @$line_aref), "\n";
        }
        else {
            $line_aref->[0] = "No_Cog";
            $line_aref->[1] = "No_Cog";
            print $COGOUT join("\t", @$line_aref), "\n";
        }
        
        #func file
        if ( $c2f_href->{$line_aref->[0]} ) {
            my $func;
            if( length($c2f_href->{$line_aref->[0]}) > 1 ) {
                my @split_func = split("",$c2f_href->{$line_aref->[0]});
                $func = $split_func[int(rand(length($c2f_href->{$line_aref->[0]})))];
            }
            else {
                $func = $c2f_href->{$line_aref->[0]};
            }
            $line_aref->[0] = $func;
            $line_aref->[1] = $func;
            print $FNOUT join("\t", @$line_aref), "\n";
        }
        else {
            $line_aref->[0] = "No_Func";
            $line_aref->[1] = "No_Func";
            print $FNOUT join("\t", @$line_aref), "\n";
        }
    }
    
    close $COGOUT;
    close $FNOUT;
    
    return 1;
}


