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

#Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'tax_file|tf=s'     => \$taxonomy_results,
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
        next if ( $key =~ qr/col_/ || $key =~ qr/cpr_/ || $key =~ qr/soil_/ );
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
        next if ( $pkey =~ qr/col_/ || $pkey =~ qr/cpr_/ || $pkey =~ qr/soil_/ );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( $ckey =~ qr/col_/ || $ckey =~ qr/cpr_/ || $ckey =~ qr/soil_/ );
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
        next if ( $pkey =~ qr/col_/ || $pkey =~ qr/cpr_/ || $pkey =~ qr/soil_/ );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( $ckey =~ qr/col_/ || $ckey =~ qr/cpr_/ || $ckey =~ qr/soil_/ );
            my $class_href = $phylum_href->{$ckey};
            foreach my $okey ( keys %$class_href ) {
                next if ( $okey =~ qr/col_/ || $okey =~ qr/cpr_/ || $okey =~ qr/soil_/ );
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
        next if ( $pkey =~ qr/col_/ || $pkey =~ qr/cpr_/ || $pkey =~ qr/soil_/ );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( $ckey =~ qr/col_/ || $ckey =~ qr/cpr_/ || $ckey =~ qr/soil_/ );
            my $class_href = $phylum_href->{$ckey};
            foreach my $okey ( keys %$class_href ) {
                next if ( $okey =~ qr/col_/ || $okey =~ qr/cpr_/ || $okey =~ qr/soil_/ );
                my $order_href = $class_href->{$okey};
                foreach my $fkey ( keys %$order_href ) {
                    next if ( $fkey =~ qr/col_/ || $fkey =~ qr/cpr_/ || $fkey =~ qr/soil_/ );
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
        next if ( $pkey =~ qr/col_/ || $pkey =~ qr/cpr_/ || $pkey =~ qr/soil_/ );
        my $phylum_href = $eu_node->{$pkey};
        foreach my $ckey ( keys %$phylum_href ) {
            next if ( $ckey =~ qr/col_/ || $ckey =~ qr/cpr_/ || $ckey =~ qr/soil_/ );
            my $class_href = $phylum_href->{$ckey};
            foreach my $okey ( keys %$class_href ) {
                next if ( $okey =~ qr/col_/ || $okey =~ qr/cpr_/ || $okey =~ qr/soil_/ );
                my $order_href = $class_href->{$okey};
                foreach my $fkey ( keys %$order_href ) {
                    next if ( $fkey =~ qr/col_/ || $fkey =~ qr/cpr_/ || $fkey =~ qr/soil_/ );
                    my $family_href = $order_href->{$fkey};
                     foreach my $gkey ( keys %$family_href ) {
                        next if ( $gkey =~ qr/col_/ || $gkey =~ qr/cpr_/ || $gkey =~ qr/soil_/ );
                        filter( $family_href->{$gkey} );
                     }
                }
            }
        }
    }
}
        

sub create_tree_node {
    my %node;
    $node{'cpr_yng'}   = 0;
    $node{'cpr_old'}   = 0;
    $node{'col_old'}   = 0;
    $node{'col_yng'}   = 0;
    $node{'soil_old'}  = 0;
    $node{'soil_old'}  = 0;
    
    return \%node;
}

sub add_part {
    my ( $part, $node ) = @_;
    if ( $part =~ qr/cpr_old/i ) {
        $node->{'cpr_old'} += 1;
    }
    elsif ( $part =~ qr/cpr_yng/i ) {
        $node->{'cpr_yng'} += 1;
    }
    elsif ( $part =~ qr/soil_old/i ) {
        $node->{'soil_old'} += 1;
    }
    elsif ( $part =~ qr/soil_yng/i ) {
        $node->{'soil_yng'} += 1;
    }
    elsif ( $part =~ qr/col_old/i ) {
        $node->{'col_old'} += 1;
    }
    elsif ( $part =~ qr/col_yng/i ) {
        $node->{'col_yng'} += 1;
    }
}

sub get_node_total {
    my ( $node ) = @_;
    my $total = $node->{'cpr_old'} + $node->{'col_old'} + $node->{'soil_old'} + $node->{'cpr_yng'} + $node->{'col_yng'} + $node->{'soil_yng'};
    return $total;
}

#filters one stage down from the href passed in
sub filter {
    my ( $href ) = @_;
    my $total = get_node_total( $href );
    foreach my $key (keys %$href) {
        next if ( $key =~ qr/col_/ || $key =~ qr/cpr_/ || $key =~ qr/soil_/ );
        if ( (get_node_total($href->{$key}))/$total < $filter ||
             get_node_total($href->{$key}) < 100 ) {
            delete $href->{$key};
        }
    }
}

sub print_line {
    my ( $href, $phynum ) = @_;
    my $col_o = $href->{'col_old'};
    my $cpr_o = $href->{'cpr_old'};
    my $soil_o = $href->{'soil_old'};
    my $col_y = $href->{'col_yng'};
    my $cpr_y = $href->{'cpr_yng'};
    my $soil_y = $href->{'soil_yng'};
    my $total = get_node_total($href);
    my $line = "\t$phynum\tcol_old=$col_o;col_yng=$col_y;cpr_old=$cpr_o;cpr_yng=$cpr_y;soil_old=$soil_o;soil_yng=$soil_y;total=$total\n";
    return $line;
}

sub print_data {
    my ( $tree_href, $FH ) = @_;
    my $eu_href = $tree_href->{'Eukaryota'};
    my $phynumber = 0;
    print $FH "K\tEukaryota", print_line( $eu_href, $phynumber );
    foreach my $pkey ( keys %$eu_href ) { #Phylum
        next if ( $pkey =~ qr/col_/ || $pkey =~ qr/cpr_/ || $pkey =~ qr/soil_/ );
        $phynumber++;
        my $phylum_href = $eu_href->{$pkey};
        print $FH "P\t$pkey", print_line( $phylum_href, $phynumber );
        foreach my $ckey ( keys %$phylum_href ) { #Class
            next if ( $ckey =~ qr/col_/ || $ckey =~ qr/cpr_/ || $ckey =~ qr/soil_/ );
            my $class_href = $phylum_href->{$ckey};
            print $FH "C\t$ckey", print_line( $class_href, $phynumber );
            foreach my $okey ( keys %$class_href ) { #Order
                next if ( $okey =~ qr/col_/ || $okey =~ qr/cpr_/ || $okey =~ qr/soil_/ );
                my $order_href = $class_href->{$okey};
                print $FH "O\t$okey", print_line( $order_href, $phynumber );
                foreach my $fkey ( keys %$order_href ) { #Family
                    next if ( $fkey =~ qr/col_/ || $fkey =~ qr/cpr_/ || $fkey =~ qr/soil_/ );
                    my $family_href = $order_href->{$fkey};
                    print $FH "F\t$fkey", print_line( $family_href, $phynumber );
                    foreach my $gkey ( keys %$family_href ) { #Genus
                        next if ( $gkey =~ qr/col_/ || $gkey =~ qr/cpr_/ || $gkey =~ qr/soil_/ );
                        my $genus_href = $family_href->{$gkey};
                        print $FH "G\t$gkey", print_line( $genus_href, $phynumber );
                        foreach my $skey ( keys %$genus_href ) { #Species
                            next if ( $skey =~ qr/col_/ || $skey =~ qr/cpr/ || $skey =~ qr/soil_/ );
                            my $species_href = $genus_href->{$skey};
                            print $FH "S\t$skey", print_line( $species_href, $phynumber );
                        }
                    }
                }
            }
        }
    } 
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
open( my $OFH, ">", $out );
print_data( \%tree, $OFH);
close( $OFH );



