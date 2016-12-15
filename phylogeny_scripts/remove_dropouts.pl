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
my $name_file;
my $gene_file;
my $drop_out_dir;
my $ass_dir;
my $out_dir; 

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'name_file=s'       => \$name_file,
            'out_dir|o=s'       => \$out_dir,
            'gene_file|gf=s'    => \$gene_file,
            'drop_dir|dd=s'      => \$drop_out_dir,
            'assem_dir=s'       => \$ass_dir,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#Subroutines
sub print_franken_genes {
    my ( $genes_w_dropouts_href, $ids_w_hash_of_full_genes_href, $out ) = @_;
    
    foreach my $gene_name ( keys %$genes_w_dropouts_href ) {
        open my $fh, ">", "$out/$gene_name.no_drops.fasta";

        foreach my $id ( keys %$ids_w_hash_of_full_genes_href ) {
            print $fh ">$id\n";
            my $id_href = $ids_w_hash_of_full_genes_href->{$id};
            my $full_gene_aref = $id_href->{$gene_name};
            my $gene_href = $genes_w_dropouts_href->{$gene_name};
            print $gene_name, "\t$id\t", scalar(@$full_gene_aref), "\t", scalar(keys %$gene_href), "\n";
            
            for ( my $i = 0; $i < scalar(@$full_gene_aref); $i++ ) { 
                if ( !exists $gene_href->{$i} ) {
                    print $fh $full_gene_aref->[$i];
                }
            }
            print $fh "\n";
        }
        close( $fh );
    }
    return 1;
}


# MAIN
my %drop_outs_per_gene;
my %hash_of_ids_with_hash_of_full_genes;

#open the drop outs per gene
my $gene_name_fo = file( $gene_file );
my @gene_names = $gene_name_fo->slurp( chomp=>1 );

my $name_file_fo = file( $name_file );
my @id_names = $name_file_fo->slurp( chomp=>1 );

#set the drop out spots for each gene
foreach my $gene ( @gene_names ) {
    my $drop_fo = file( "$drop_out_dir/$gene.dropout_pos.txt" );
    my @lines = $drop_fo->slurp( chomp=>1 );
    $drop_outs_per_gene{$gene} = {};
    my $href = $drop_outs_per_gene{$gene};
    
    foreach my $line ( @lines ) {
        if ( $line =~ qr/[0-9]+,/ ) {
            my @drops = split qr/,/, $line;
            foreach my $drop_out (@drops) {
                if ( !exists $href->{$drop_out} ) {
                    $href->{$drop_out} = 0;
                }
            }
        }
    }
}

foreach my $id ( @id_names ) {
    my $id_ass_fo = file( "$ass_dir/$id.asm_genes.fasta" );
    my @data_lines = $id_ass_fo->slurp( chomp=>1 );
    my $gn; #Gene Name
    $hash_of_ids_with_hash_of_full_genes{$id} = {};
    my $href = $hash_of_ids_with_hash_of_full_genes{$id};
    
    foreach my $line ( @data_lines ) {
        if ( $line =~ qr/^>/ ) {
            $gn = $line;
            $gn =~ s/>//;
            $href->{$gn} = [];
        }
        else {
            #add nucleotides to full gene
            if ( $line =~ qr/\A[actg]+\z/i ) {
                my @nucleotides = split qr//, $line;
                foreach my $snclt ( @nucleotides ) {
                    push @{$href->{$gn}}, $snclt;
                }
            }
        }
    }
}
print_franken_genes( \%drop_outs_per_gene, \%hash_of_ids_with_hash_of_full_genes, $out_dir);






