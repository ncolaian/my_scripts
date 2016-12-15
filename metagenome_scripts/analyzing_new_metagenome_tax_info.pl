#! usr/bin/env perl

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

# Variables
my $help = 0;
my $man = 0;
my $taxonomy_results;
my $out;
my $filter;
my $sample_ids;

#Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'tax_file|tf=s'     => \$taxonomy_results,
            'sample_ids|si=s'   => \$sample_ids,
            'out|o=s'           => \$out,
            'read_filt|rf:i'    => \$filter,   
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Handle filter that is not specified and when it is specified
# This param just makes sure the stuff we pick is significant
if ( !defined $filter ) {
    $filter = .01;
}
else {
    $filter = $filter/100;
}

#hash that represents the tree variable
my %tree;

#sample_id_array
my $sample_fo = file( $sample_ids );
my @sample_ids = $sample_fo->slurp( chomp=>1 );

#
my %check_hash;
foreach my $id ( @sample_ids ) {
    $check_hash{$id} = 1;
}

# Subroutines

sub count_eurkaryotes {
    my ( $arefs_tax_line_aref, $tree_href ) = @_;
    $tree_href->{'Eukaryota'} = create_tree_node();
    my $node_href = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar(@$arefs_tax_line_aref); $i++ ) {
        my $line_aref = $arefs_tax_line_aref->[$i];
        if ( $line_aref->[1] =~ qr/eukaryota/i ) {
            add_part( $line_aref->[0] ,$node_href);
        }
    }
}

sub count_and_filter_phylum {
    my ( $arefs_tax_line_aref, $part_of_tax, $tree_href ) = @_;
    my $eu_node = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar( @$arefs_tax_line_aref ); $i++ ) {
        next if ( $arefs_tax_line_aref->[$i] == 1 );
        my $line_aref = $arefs_tax_line_aref->[$i];
        my @split_tax = split /;/, $line_aref->[1];
        if ( scalar ( @split_tax ) == 1 ) { #Gets rid of lines that are fully defined
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        
        if ( !defined $eu_node->{ $split_tax[$part_of_tax] } ) {
            $eu_node->{ $split_tax[$part_of_tax] } = create_tree_node();
        }
        add_part( $line_aref->[0], $eu_node->{ $split_tax[$part_of_tax] } );
    }
    filter( $eu_node );
}

sub count_and_filter_class {
    my ( $arefs_tax_line_aref, $part_of_tax, $tree_href ) = @_;
    my $eu_node = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar( @$arefs_tax_line_aref ); $i++ ) {
        next if ( $arefs_tax_line_aref->[$i] == 1 );
        my $line_aref = $arefs_tax_line_aref->[$i];
        my @split_tax = split /;/, $line_aref->[1];
        if ( scalar ( @split_tax ) == 2 ) { #Gets rid of lines that are fully defined
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        #Handles getting rid of phylum lines that we don't need to worry about
        if ( !defined $eu_node->{ $split_tax[$part_of_tax-1] } ||
             scalar( @split_tax ) <= $part_of_tax ) {
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        else {
            my $phylum_href = $eu_node->{ $split_tax[$part_of_tax-1] };
            if ( !defined $phylum_href->{ $split_tax[$part_of_tax] } ) {
                $phylum_href->{ $split_tax[$part_of_tax] } = create_tree_node();
            }
            add_part( $line_aref->[0], $phylum_href->{ $split_tax[$part_of_tax] } );
        }
    }
    foreach my $key ( keys %$eu_node ) {
        next if ( exists $check_hash{$key} );
        filter( $eu_node->{$key} );
    }
}

sub count_and_filter_order {
    my ( $arefs_tax_line_aref, $part_of_tax, $tree_href ) = @_;
    my $eu_node = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar( @$arefs_tax_line_aref ); $i++ ) {
        next if ( $arefs_tax_line_aref->[$i] == 1 );
        my $line_aref = $arefs_tax_line_aref->[$i];
        my @split_tax = split /;/, $line_aref->[1];
        if ( scalar ( @split_tax ) == 3 ) { #Gets rid of lines that are fully defined
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        my $phylum_href = $eu_node->{ $split_tax[$part_of_tax-2] };
        #Handles getting rid of class lines that we don't need to worry about
        if ( !defined $phylum_href->{ $split_tax[$part_of_tax-1] } ||
             scalar( @split_tax ) <= $part_of_tax ) {
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        else {
            my $class_href = $phylum_href->{ $split_tax[$part_of_tax-1] };
            if ( !defined $class_href->{ $split_tax[$part_of_tax] } ) {
                $class_href->{ $split_tax[$part_of_tax] } = create_tree_node();
            }
            add_part( $line_aref->[0], $class_href->{ $split_tax[$part_of_tax]} );
        }
    }
    foreach my $pkey ( keys %$eu_node ) {
        next if ( exists $check_hash{$pkey} );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( exists $check_hash{$ckey} );
            filter( $phylum_href->{$ckey} );
        }
    }
}

sub count_and_filter_family {
    my ( $arefs_tax_line_aref, $part_of_tax, $tree_href ) = @_;
    my $eu_node = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar( @$arefs_tax_line_aref ); $i++ ) {
        next if ( $arefs_tax_line_aref->[$i] == 1 );
        my $line_aref = $arefs_tax_line_aref->[$i];
        my @split_tax = split /;/, $line_aref->[1];
        if ( scalar ( @split_tax ) == 4 ) { #Gets rid of lines that are fully defined
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        my $phylum_href = $eu_node->{ $split_tax[$part_of_tax-3] };
        my $class_href = $phylum_href->{ $split_tax[$part_of_tax-2] };
        if ( !defined $class_href->{ $split_tax[$part_of_tax-1] } ||
             scalar( @split_tax ) <= $part_of_tax ) {
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        else {
            my $order_href = $class_href->{ $split_tax[$part_of_tax-1] };
            if ( !defined $order_href->{ $split_tax[$part_of_tax] } ) {
                $order_href->{ $split_tax[$part_of_tax] } = create_tree_node();
            }
            add_part( $line_aref->[0], $order_href->{ $split_tax[$part_of_tax]} );
        }
    }
    foreach my $pkey ( keys %$eu_node ) {
        next if ( exists $check_hash{$pkey} );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( exists $check_hash{$ckey} );
            my $class_href = $phylum_href->{$ckey};
            foreach my $okey ( keys %$class_href ) {
                next if ( exists $check_hash{$okey} );
                filter( $class_href->{$okey} );
            }
        }
    }
}

sub count_and_filter_genus {
    my ( $arefs_tax_line_aref, $part_of_tax, $tree_href ) = @_;
    my $eu_node = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar( @$arefs_tax_line_aref ); $i++ ) {
        next if ( $arefs_tax_line_aref->[$i] == 1 );
        my $line_aref = $arefs_tax_line_aref->[$i];
        my @split_tax = split /;/, $line_aref->[1];
        if ( scalar ( @split_tax ) == 5 ) { #Gets rid of lines that are fully defined
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        my $phylum_href = $eu_node->{ $split_tax[$part_of_tax-4] };
        my $class_href = $phylum_href->{ $split_tax[$part_of_tax-3] };
        my $order_href = $class_href->{ $split_tax[$part_of_tax-2] };
        if ( !defined $order_href->{ $split_tax[$part_of_tax-1] } ||
             scalar( @split_tax ) <= $part_of_tax ) {
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        else {
            my $family_href = $order_href->{ $split_tax[$part_of_tax-1] };
            if ( !defined $family_href->{ $split_tax[$part_of_tax] } ) {
                $family_href->{ $split_tax[$part_of_tax] } = create_tree_node();
            }
            add_part( $line_aref->[0], $family_href->{ $split_tax[$part_of_tax]} );
        }
    }
    foreach my $pkey ( keys %$eu_node ) {
        next if ( exists $check_hash{$pkey} );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( exists $check_hash{$ckey} );
            my $class_href = $phylum_href->{$ckey};
            foreach my $okey ( keys %$class_href ) {
                next if ( exists $check_hash{$okey} );
                my $order_href = $class_href->{$okey};
                foreach my $fkey ( keys %$order_href ) {
                    next if ( exists $check_hash{$fkey} );
                    filter( $order_href->{$fkey} );
                }
            }
        }
    }
}

sub count_and_filter_species {
    my ( $arefs_tax_line_aref, $part_of_tax, $tree_href ) = @_;
    my $eu_node = $tree_href->{'Eukaryota'};
    for ( my $i = 1; $i < scalar( @$arefs_tax_line_aref ); $i++ ) {
        next if ( $arefs_tax_line_aref->[$i] == 1 );
        my $line_aref = $arefs_tax_line_aref->[$i];
        my @split_tax = split /;/, $line_aref->[1];
        if ( scalar ( @split_tax ) == 6 ) { #Gets rid of lines that are fully defined
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        
        my $phylum_href = $eu_node->{ $split_tax[$part_of_tax-5] };
        my $class_href = $phylum_href->{ $split_tax[$part_of_tax-4] };
        my $order_href = $class_href->{ $split_tax[$part_of_tax-3] };
        my $family_href = $order_href->{ $split_tax[$part_of_tax-2] };
        if ( !defined $family_href->{ $split_tax[$part_of_tax-1] } ||
             scalar( @split_tax ) <= $part_of_tax ) {
            $arefs_tax_line_aref->[$i] = 1;
            next;
        }
        else {
            my $genus_href = $family_href->{ $split_tax[$part_of_tax-1] };
            if ( !defined $genus_href->{ $split_tax[$part_of_tax] } ) {
                $genus_href->{ $split_tax[$part_of_tax] } = create_tree_node();
            }
            add_part( $line_aref->[0], $genus_href->{ $split_tax[$part_of_tax]} );
        }
    }
    foreach my $pkey ( keys %$eu_node ) {
        next if ( exists $check_hash{$pkey} );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( exists $check_hash{$ckey} );
            my $class_href = $phylum_href->{$ckey};
            foreach my $okey ( keys %$class_href ) {
                next if ( exists $check_hash{$okey} );
                my $order_href = $class_href->{$okey};
                foreach my $fkey ( keys %$order_href ) {
                    next if ( exists $check_hash{$fkey} );
                    my $family_href = $order_href->{$fkey};
                     foreach my $gkey ( keys %$family_href ) {
                        next if ( exists $check_hash{$gkey} );
                        filter( $family_href->{$gkey} );
                     }
                }
            }
        }
    }
}
        

sub create_tree_node {
    my %node;
    foreach my $id ( @sample_ids ) {
        $node{$id} = 0;
    }
    
    return \%node;
}

sub add_part {
    my ( $part, $node ) = @_;
    $node->{$part} += 1;
}

sub get_node_total {
    my ( $node ) = @_;
    my $total = 0;
    foreach my $id ( @sample_ids ) {
        $total += $node->{$id};
    }
    return $total;
}

#filters one stage down from the href passed in
sub filter {
    my ( $href ) = @_;
    my $total = get_node_total( $href );
    foreach my $key (keys %$href) {
        next if ( exists $check_hash{$key} );
        if ( (get_node_total($href->{$key}))/$total < $filter ||
             get_node_total($href->{$key}) < 100 ) {
            delete $href->{$key};
        }
    }
}

sub print_line {
    my ( $href ) = @_;
    next if ( ref( $href ) ne "HASH");
    my $total = get_node_total($href);
    my $line;
    foreach my $id ( @sample_ids ) {
        $line .= $href->{$id};
        $line .=",";
    }
    $line .= "$total\n";
    return $line;
}

sub print_data {
    my ( $tree_href, $FH ) = @_;
    my $eu_href = $tree_href->{'Eukaryota'};
    print $FH get_id_header();
    print $FH "K,Eukaryota,NA,", print_line( $eu_href);
    foreach my $pkey ( keys %$eu_href ) { #Phylum
        next if ( exists $check_hash{$pkey} );
        my $phylum_href = $eu_href->{$pkey};
        print $FH "P,$pkey,Eukaryota,", print_line( $phylum_href );
        foreach my $ckey ( keys %$phylum_href ) { #Class
            next if ( exists $check_hash{$ckey} );
            my $class_href = $phylum_href->{$ckey};
            print $FH "C,$ckey,$pkey,", print_line( $class_href );
            foreach my $okey ( keys %$class_href ) { #Order
                next if ( exists $check_hash{$okey} );
                my $order_href = $class_href->{$okey};
                print $FH "O,$okey,$ckey,", print_line( $order_href );
                foreach my $fkey ( keys %$order_href ) { #Family
                    next if ( exists $check_hash{$fkey} );
                    my $family_href = $order_href->{$fkey};
                    print $FH "F,$fkey,$okey,", print_line( $family_href );
                    foreach my $gkey ( keys %$family_href ) { #Genus
                        next if ( exists $check_hash{$gkey} );
                        my $genus_href = $family_href->{$gkey};
                        print $FH "G,$gkey,$fkey,", print_line( $genus_href );
                        foreach my $skey ( keys %$genus_href ) { #Species
                            next if ( $skey =~ qr/col_/ || $skey =~ qr/cpr/ || $skey =~ qr/soil_/ );
                            my $species_href = $genus_href->{$skey};
                            print $FH "S,$skey,$gkey,", print_line( $species_href );
                        }
                    }
                }
            }
        }
    } 
}

sub get_id_header {
    my $header = "Taxa_level,Group_Name,Group_Above,";
    $header .= join ",", @sample_ids;
    $header .=",Total\n";
    return $header;
}

# Main
my $tax_fo = file( $taxonomy_results );
my @arefs_tax_slurp = $tax_fo->slurp( chomp=>1, split=>qr/\t/ );

#count eukaryotes and place a 1 where the line is removed.
count_eurkaryotes( \@arefs_tax_slurp, \%tree);

#count phylum, and then filter them, and add to tree
# 1 = phylum, 2 = class, 3 = order, 4 = family, 5 = genus, 6 = species
count_and_filter_phylum( \@arefs_tax_slurp, 1, \%tree );
count_and_filter_class(  \@arefs_tax_slurp, 2, \%tree );
count_and_filter_order(  \@arefs_tax_slurp, 3, \%tree );
count_and_filter_family( \@arefs_tax_slurp, 4, \%tree );
count_and_filter_genus(  \@arefs_tax_slurp, 5, \%tree );
count_and_filter_species(\@arefs_tax_slurp, 6, \%tree );
open ( my $OFH, ">", $out );
print_data( \%tree, $OFH);
close( $OFH );



