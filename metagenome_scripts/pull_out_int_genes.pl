#! /usr/bin/env perl

# This program will take the output from the local r script find_adj_val.R, which reports the cogs that have adjacent genomes that have different abundance calls. It will then create a directory containing fasta files of all the genes that are associated with that cog in the genome.

#Need to have a way to 

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;
use Readonly;

#Create the hash that will be used to translate DNA to Protein Sequence
Readonly::Hash my %CODON_TBL => ('TCA'=>'S','TCC'=>'S','TCG'=>'S',
        'TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y',
        'TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_',
        'TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P',
        'CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q',
        'CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I',
        'ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T',
        'ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S',
        'AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V',
        'GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D',
     'GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G',
        'GGT'=>'G');

#My Variables
my $help = 0;
my $man = 0;
my $dafe_dir;
my $out_dir;
my $adj_diff_file;
my $ordered_genome_file;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'dafe_dir|dd=s'     =>  \$dafe_dir,
            'out_dir|o=s'       =>  \$out_dir,
            'genome_file|gf=s'  =>  \$ordered_genome_file,
            'diff_file|diff=s'  =>  \$adj_diff_file,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##

$logger->info("Reading in the files, and chhecking to make sure the DAFE directory exists");

if ( !-d $dafe_dir) { croak "DAFE directory passed does not exist" }

my $diff_fo = file($adj_diff_file);
my @diff_file = $diff_fo->slurp( split=>qr/\t/, chomp=>1 );

my $ordered_fo = file( $ordered_genome_file );
my @ord_genome = $ordered_fo->slurp( chomp=>1 );


$logger->info( "Creating a hash containing the genomes to pull the cog gene sequences from" );

my $hash_of_cogs_to_genomes_href = create_cogs_to_genomes_href( \@diff_file, \@ord_genome );

$logger->info( "Creating Directory Structure and Fastas" );

make_dir_structure( $hash_of_cogs_to_genomes_href );

create_fastas( $hash_of_cogs_to_genomes_href );

perform_muscle_alignment();



## Subroutines ##

sub create_cogs_to_genomes_href {
    my ( $diff_aref, $ordered_aref ) = @_;
    
    #Href containing the cog to genomes - genomes are separated by a tab
    my %href;
    my $throw_header = shift @$diff_aref;
    
    foreach my $line_aref ( @$diff_aref ) {
        my $cog = shift @$line_aref;
        my $up_g = $line_aref->[0];
        my $genomes = join( "\t", @$line_aref );
        #get the outgroup that is 20 genomes away from this
        my $out = get_outgroup($up_g, $ordered_aref);
        my $g = join("\t", @$out);
        $genomes .= "\t$g";
        
        #set or add to cog
        if ( $href{$cog} ) {
            $href{$cog} = $href{$cog} . "\t$genomes";
        }
        else {
            $href{$cog} = $genomes;
        }
    }
    
    return \%href;
}

sub get_outgroup {
    my ( $genome, $ord_aref ) = @_;
    my @out_genome;
    for ( my $i = 0; $i < scalar(@$ord_aref); $i++ ) {
        if ( $ord_aref->[$i] eq $genome ) {
            #Get all the genomes around this interesting pair
            push @out_genome, $ord_aref->[$i+21];
            push @out_genome, $ord_aref->[$i-20];
            push @out_genome, $ord_aref->[$i-1];
            push @out_genome, $ord_aref->[$i-2];
            push @out_genome, $ord_aref->[$i-3];
            push @out_genome, $ord_aref->[$i+2];
            push @out_genome, $ord_aref->[$i+3];
            push @out_genome, $ord_aref->[$i+4];
            last;
        }
    }
    if ( !@out_genome ) { croak "$genome does not have an out genome"}
    
    return \@out_genome;
}

sub make_dir_structure {
    my ( $cog_href ) = @_;
    if ( !-d $out_dir ) {
        mkdir $out_dir;
    }
    
    foreach my $keys ( keys %$cog_href ) {
        mkdir "$out_dir/$keys";
        mkdir "$out_dir/$keys";
        mkdir "$out_dir/$keys";
    }
    
    return 1;
}

sub create_fastas {
    my ( $cog_href ) = @_;
    
    foreach my $cog ( keys %$cog_href ) {
        my @genomes = split /\t/, $cog_href->{$cog};
        
        #open fasta file to be written
        open my $OUT, ">", "$out_dir/$cog/pre_nucl_$cog.fasta";
        
        foreach my $g ( @genomes ) {
            my $gene_ids_href = get_genes( $g, $cog);
            
            #read in fasta file
            my $fasta_fo = file("$dafe_dir/$g/$g.fna");
            my @fasta = $fasta_fo->slurp( chomp=>1 );
            
            #write fasta files
            my $write = 0; #holds if we want the sequence or not after name
            my @sequence; #holds current sequence
            my $st_en_aref; #holds the start and stop positions of the genes
            my $prev_line;
            foreach my $line ( @fasta ) {
                #deals with name lines
                if ( $line =~ qr/>/ ) {
                    $line =~ s/>//;
                    if ( $write == 1 ) {
                        $write = 0;
                        my $seq = join "", @sequence;
                        my @seq_arr = split //, $seq;
                        #handles multiple genes on a single contig
                        my $num = 1;
                        foreach my $pos ( @$st_en_aref) {
                            my ($start, $end) = split /\t/, $pos;
                            print $OUT ">$g\_$prev_line\_$num\n";
                            for ( my $i = ($start-1); $i < $end; $i++ ) {
                                print $OUT $seq_arr[$i];
                            }
                            print $OUT "\n";
                        }
                        undef $st_en_aref;
                        undef $prev_line;
                    }
                    if ( $gene_ids_href->{$line} ) {
                        $prev_line = $line;
                        $write = 1;
                        $st_en_aref = $gene_ids_href->{$line};
                    }
                }
                #deals with sequence lines
                else {
                    if ( $write == 1 ) {
                        push @sequence, $line;
                    }
                }
            }
            #This deals with the last line
            if ( $write == 1 ) {
                $write = 0;
                my $seq = join "", @sequence;
                my @seq_arr = split //, $seq;
                #handles multiple genes on a single contig
                my $num = 1;
                foreach my $pos ( @$st_en_aref) {
                    my ($start, $end) = split /\t/, $pos;
                    print $OUT ">$g\_$prev_line\_$num\n";
                    for ( my $i = ($start-1); $i < $end; $i++ ) {
                        print $OUT $seq_arr[$i];
                    }
                    print $OUT "\n";
                }
            }        
        }
        close $OUT;
        
        #Need to make sure the fasta sequences are all on the same strand - Using my pre_mult_align_seq_flip.pl
        my $cmd = "perl /nas02/home/n/c/ncolaian/my_scripts/utility_scripts/pre_mult_align_seq_flip.pl -af $out_dir/$cog/pre_nucl_$cog.fasta -o $out_dir/$cog/nucl_$cog.fasta";
        $logger->debug( $cmd );
        system( $cmd );
        
        #Need to create a protein file
        translate_file( "$out_dir/$cog/nucl_$cog.fasta", $cog );
    }
    return 1;
}

sub get_genes {
    my ( $genome, $c ) = @_;
    
    # First need to get ID from all_annote file
    my $annote_fo = file("$dafe_dir/$genome/all_annote.txt");
    my @annote = $annote_fo->slurp( chomp=>1, split=>qr/\t/ );
    
    #Go through and get ID's for all the things that equal cogs
    my $header_aref = shift @annote;
    my $cog_col = get_cog_col( $header_aref ); #get cog col number
    
    my %ids;
    foreach my $line_aref ( @annote ) {
        if ( $line_aref->[$cog_col] =~ qr/$c/i ) {
            $ids{$line_aref->[0]} = 1;
        }
    }
    
    #Now need to get gene names from gff file - link to sequences in fna file
    my $gff_fo = file("$dafe_dir/$genome/$genome.gff");
    my @gff = $gff_fo->slurp( chomp=>1, split=>qr/\t/ );
    
    my %gene_ids;
    foreach my $aref ( @gff ) {
        my $line_with_id = $aref->[scalar(@$aref)-1];
        if ( $line_with_id !~ qr/ID/ ) {
            next;
        }

        if ( contain_an_id( \%ids, $line_with_id) ) {
            my $start = $aref->[3];
            my $end = $aref->[4];
            #takes into account multiple genes in a contig
            if ( $gene_ids{$aref} ) {
                $gene_ids{$aref->[0]} = push @{$gene_ids{$aref->[0]}}, ["$start\t$end"];
            }
            else {
                $gene_ids{$aref->[0]} = ["$start\t$end"];
            }
        }
    }
    
    return \%gene_ids;
}

sub get_cog_col {
    my ( $line_aref ) = @_;
    my $col_num = 0;
    foreach my $col ( @$line_aref ) {
        if ( $col =~ m/cog/ ) {
            last;
        }
        $col_num++;
    }
    return $col_num;
}

sub contain_an_id {
    my ( $href, $line ) = @_;
    #print Dumper($href), "hi\n";
    
    foreach my $key ( keys %$href ) {
        #print "$line\t$key\n";
        if ( $line =~ qr/$key/i ) {
            return 1;
        }
    }
    return 0;
}

sub translate_file {
    my ( $file_path, $cur_cog ) = @_;
    open my $IN, "<", $file_path;
    open my $OUT, ">", "$out_dir/$cur_cog/prot_$cur_cog.fasta";
    
    my $prev_gene;
    foreach my $line ( <$IN> ) {
        if ( $line =~ /^>/ ) {
            print $OUT $line;
            $prev_gene = $line;
            next;
        }
        chomp $line;
        
        #Check to make sure string is divisible by three. If not raise warning
        if ( (length($line) % 3) != 0 ) {
            carp "$prev_gene is not devisible by three and the trailing bases will be ignored";
        }
        
        #Check to make sure the fist codon is a start
        my $first_codon = substr $line, 0, 3;
        if( $CODON_TBL{$first_codon} ne "M" ) {
            carp "$prev_gene has been reverse compimented becasue the gene does not start with a Methionine";
            $line = reverse_comp( $line );
        }
        
        #Translate DNA to Prot
        my $transline= "";
        for( my $i = 0; $i < (length($line)-2); $i = $i + 3 ) {
            my $codon = substr $line, $i, 3;
            $codon = uc $codon;
            
            if ( $CODON_TBL{$codon} ) {
                $transline .= $CODON_TBL{$codon};
            }
            elsif ( $codon =~ m/N/g ) {
                $transline .= 'X';
            }
            else {
                croak "There is an unknown codon in $prev_gene";
            }
        }
        
        # remove any trailing stop codon symbols (ie "_")
        $transline =~ s/_$//;
        
        print $OUT "$transline\n";
    }
    close $IN;
    close $OUT;
}

sub reverse_comp {
    my ( $str ) = @_;
    $str = reverse $str;
    $str =~ tr/ACGTacgt/TGCAtgca/;
    return $str;
}

sub perform_muscle_alignment {
    opendir(my $DIR, $out_dir);
    
    my @cogs;
    foreach my $dir ( readdir $DIR ) {
        if ( !-d $dir ) {
            next;
        }
        my $prot_cmd = "bsub -o $out_dir/muscle.out muscle -in $out_dir/$dir/prot_$dir.fasta -out prot_align_$dir.fasta";
        my $nucl_cmd = "bsub -o $out_dir/muscle.out muscle -in $out_dir/$dir/nucl_$dir.fasta -out nucl_align_$dir.fasta";
        $logger->debug( $prot_cmd );
        $logger->debug( $nucl_cmd );
        system( $prot_cmd );
        system( $nucl_cmd );
    }
}


__END__
=head1 Pull_out_int_genes.pl

This program will pull out the genes from associated with cogs from interesting genome pairs.

=head1 VERSION

0.1

=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT

    

=head1 DEPENDENCIES

    
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2017, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut
