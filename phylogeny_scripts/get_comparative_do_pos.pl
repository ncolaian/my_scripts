#! usr/bin/evn perl

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
my $name_file;
my $gene_file; 
my $out_dir;
my $ref_gen_file;
my $drop_out_dir;
my $ass_dir;


# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'name_file=s'       => \$name_file,
            'out_dir|o=s'       => \$out_dir,
            'gene_name|gn=s'    => \$gene_file,
            'drop_dir|dd=s'     => \$drop_out_dir,
            'ref_file|ref=s'    => \$ref_gen_file,
            'assem_dir=s'       => \$ass_dir,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Subroutines that will occur in the main program

sub create_alignment {
    my ( $ass_dir, $ref_gene_file, $out_dir, $name_aref ) = @_;
    foreach my $id (@$name_aref) {
        my $id_file_ass = "$ass_dir/$id.asm_genes.fasta";
        my $cmd = "bsub -o $out_dir/lsf.out /nas02/home/j/r/jrwang/local/bin/blasr $id_file_ass $ref_gene_file -sam -out $out_dir/$id.sam";
        system($cmd);
    }
    return 1;
}

sub wait_for_bsubs {
    my $tempdir = tempdir();
    my ($fh, $filename) = tempfile();
    my $bjobs = `bjobs -w`;
    print $fh $bjobs;
    close $fh;
    
    my $tfo = file($filename);
    my @lines = $tfo->slurp(chomp=>1);
    
    if (grep qr/blasr/, @lines) {
        sleep 60;
        return 0;
    }
    return 1;
}

sub get_gene_starts {
    my ( $out_dir, $name_aref, $gene_aref, $viril_st_href ) = @_;
    foreach my $gene ( @$gene_aref ) {
        $viril_st_href->{$gene} = {};
    }
    foreach my $id ( @$name_aref ) {
        my $file_obj = file("$out_dir/$id.sam");
        my @sam_lines = $file_obj->slurp(chomp=>1);
        foreach my $line (@sam_lines) {
            next if ($line =~ qr/@/ || $line =~ qr/^$/);
            my @align = split qr/\t/,$line;
            my $href = $viril_st_href->{$align[2]};
            $href->{$id} = $align[3];
        }
    }
}

sub find_max_adjustment_forid_in_genes {
    my ($starts_href) = @_;
    my %adjustments;
    
    foreach my $gene ( keys %$starts_href ) {
        my $href = $starts_href->{$gene};
        if ( scalar(keys %$href) < 27 ) {
            print scalar( keys %$href ), "\n";
            print "$gene has less than 27 genomes represented and should be excluded\n";
            delete $starts_href->{$gene};
            next;
        }
        my @array;
        foreach my $id ( keys %$href ) {
            push @array, $href->{$id};
        }
        my $max = max @array;
        $adjustments{$gene} = $max;
    }
    return \%adjustments;
}

sub get_adjusted_dropout_pos {
    my ( $drops_per_gene_href, $name_aref, $starts_href ) = @_;
    my %conserved_drops;
    foreach my $gene (keys %$drops_per_gene_href) {
        if ( !exists $starts_href->{$gene} ) {
            next;
        }
        my $href = $drops_per_gene_href->{$gene};
        my $s_href = $starts_href->{$gene};
        $conserved_drops{$gene} = {};
        my $drop_href = $conserved_drops{$gene};
        
        for (my $name_i = 0; $name_i < scalar( @$name_aref ); $name_i++) {
            my $drop_aref = $href->{$name_aref->[$name_i]};
            my $adjust = $s_href->{$name_aref->[$name_i]};
            
            for (my $i = 0; $i < scalar( @$drop_aref ); $i++) {
                $drop_aref->[$i] += ($adjust-1); #b/c drop_outs start at 0 position not 1 like in sam format
                if ( !exists $drop_href->{$drop_aref->[$i]} ) {
                    $drop_href->{$drop_aref->[$i]} = 0;
                }
            }
        }
    }
    return \%conserved_drops #hashref of a hashref containing dropout pos
}

sub print_conserved_drop_out_positions {
    my ( $drop_href, $out ) = @_;
    open my $fh, ">", "$out/viril_aligned_drop_sum.txt";
    foreach my $gene ( %$drop_href ) {
        next if ( !exists $drop_href->{$gene} );
        my $href = $drop_href->{$gene};
        print $fh $gene, "\tNumber of Drops = ", scalar( keys %$href), "\n";
        my @array;
        foreach my $drops ( sort { $a <=> $b} keys %$href) {
            push @array, $drops;
        }
        print $fh join(", ", @array), "\n";
    }
    close($fh);
}

sub get_full_genes {
    my ( $ids_aref, $ass_dir ) = @_;
    my %hash_of_ids_with_hash_of_full_genes;
    
    foreach my $id ( @$ids_aref ) {
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
    return \%hash_of_ids_with_hash_of_full_genes;
}

sub print_drop_out_removed_genes {
    my ( $out, $full_genes_href, $drop_href, $max_href, $adjust_href ) = @_;
    foreach my $gene ( %$max_href ) {
        next if ( !exists $drop_href->{$gene} );
        #mkdir ("$out/genes_wo_dropouts/")
        open my $fh, ">", "$out/genes_wo_dropouts/$gene.align_no_drops.fa";
        my $href = $drop_href->{$gene};
        my $href_adjust = $adjust_href->{$gene};
        my $max_adjustment = $max_href->{$gene};
        my %no_drop_genes;
        my @lengths;
        
        foreach my $id ( keys %$full_genes_href ) {
            my @gene;
            my $full_genes_id_href = $full_genes_href->{$id};
            my $full_gene_aref = $full_genes_id_href->{$gene};
            my $id_adjust = $href_adjust->{$id};
            
            for (my $i = 0; $i < scalar(@$full_gene_aref); $i++) {
                if ( !exists $href->{ ($i + ($id_adjust - 1) ) }
                    && ($i + ($id_adjust-1)) >= ($max_adjustment-1) ) {
                    push @gene, $full_gene_aref->[$i];
                    my $add = $i+ ($id_adjust-1);
                }
            }
            $no_drop_genes{$id} = \@gene;
            push @lengths, scalar(@gene);
        }
        my $min_length = min(@lengths);
        foreach my $keys ( keys %no_drop_genes ) {
            print $fh ">$keys\t$min_length\n";
            my $gene = $no_drop_genes{$keys};
            for (my $i = 0; $i < $min_length; $i++) {
                print $fh $gene->[$i];
                
            }
            print $fh "\n";
        }
        close($fh);
    }
    return 1;
}


# MAIN
my %drop_outs_per_gene; # A hash containing a hashref of genes containing ids with an aref of dropouts.
my %viril_starts; # genes with arefs of starts of the ids in the order of the gene file

#open the drop outs per gene
my $gene_name_fo = file( $gene_file );
my @gene_names = $gene_name_fo->slurp( chomp=>1 );

my $name_file_fo = file( $name_file );
my @id_names = $name_file_fo->slurp( chomp=>1 );

# Holds the dropouts
foreach my $gene ( @gene_names ) {
    my $drop_fo = file( "$drop_out_dir/$gene.dropout_pos.txt" );
    my @lines = $drop_fo->slurp( chomp=>1 );
    $drop_outs_per_gene{$gene} = {};
    my $href = $drop_outs_per_gene{$gene};
    my $id;
    
    foreach my $line ( @lines ) {
        next if $line =~ qr/^$/;
        if ( $line !~ qr/[0-9]+,/ ) {
            my @name_line_split = split qr/\t/, $line;
            $href->{$name_line_split[0]} = [];
            $id = $name_line_split[0];
        }
        else {
            my @drops = split qr/,/, $line;
            foreach my $drop_out (@drops) {
                push @{$href->{$id}}, $drop_out;
            }
        }
    }
}
create_alignment( $ass_dir, $ref_gen_file, $out_dir, \@id_names );
if (wait_for_bsubs() eq '0') {
    wait_for_bsubs();
}
get_gene_starts( $out_dir, \@id_names, \@gene_names, \%viril_starts);
my $max_values_href = find_max_adjustment_forid_in_genes( \%viril_starts );
my $conserved_dropouts_href = get_adjusted_dropout_pos( \%drop_outs_per_gene, \@id_names, \%viril_starts);
print_conserved_drop_out_positions( $conserved_dropouts_href, $out_dir );
# Get hash of ids that hold its genes with an array of the sequence
my $full_genes_href = get_full_genes( \@id_names, $ass_dir );
print_drop_out_removed_genes( $out_dir, $full_genes_href, $conserved_dropouts_href, $max_values_href, \%viril_starts );







