#! /usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# My Variables
my $help = 0;
my $man = 0;
my $dafe_dir;
my $ordered_genomes;
my $ordered_kogs;
my $group_name;
my $out_dir;

# Read in the variable from the command line
GetOptions ( 'man'  =>  \$man,
            'help'  =>  \$help,
            'ord_kogs|ok=s'     =>  \$ordered_kogs,
            'ord_genomes|go=s'  =>  \$ordered_genomes,
            'group|g=s'         =>  \$group_name,
            'dafe_dir|dd=s'     =>  \$dafe_dir,
            'out_dir|od=s'      =>  \$out_dir,
            ) || die("There was an error in the command line arguements \n");

# Pod Usage for the manual and help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

## MAIN ##
$logger->info("Create array of ordered genomes and kogs");
my $k_file = file($ordered_kogs);
my $g_file = file($ordered_genomes);

my @ord_kogs = $k_file->slurp(chomp=>1);
my @ord_genomes = $g_file->slurp(chomp=>1);

# hash that will hold hashes of kog counts per genome
$logger->info("Creating the kog count hash");
my %hash_of_hashrefs_containing_kog_counts;
foreach my $genome ( @ord_genomes ) {
    $hash_of_hashrefs_containing_kog_counts{$genome} = get_kog_count_hash($genome);
}

print_out_kog_count_table(\%hash_of_hashrefs_containing_kog_counts,\@ord_genomes, \@ord_kogs);

# R code to create the heat map
#create_heatmap($r_source);




## SUBS ##

sub get_kog_count_hash {
    my ( $genome ) = @_;
    open my $FH, "<", "$dafe_dir/$genome/all_annote.txt";
    #find the kog column
    my $kog_column = 0;
    my $first_line = <$FH>;
    my @split_fl = split /\t/, $first_line;
    foreach my $line ( @split_fl ) {
        if ( $line =~ qr/^$group_name$/ ) {
            last;
        }
        $kog_column++;
    }
    
    #go through file and count KOG info
    my %kog_counts;
    foreach my $ann_line ( <$FH> ) {
        chomp $ann_line;
        my @split = split /\t/, $ann_line;
        next if ( $split[$kog_column] =~ qr/NA/i );
        
        #check if kog existed if it does add to the current count. If not initiate
        if ( $kog_counts{$split[$kog_column]} ) {
            $kog_counts{$split[$kog_column]} = $kog_counts{$split[$kog_column]} + 1;
        }
        else {
            $kog_counts{$split[$kog_column]} = 1;
        }
    }
    
    close($FH);
    return \%kog_counts;
}

sub print_out_kog_count_table {
    my ( $count_href_hrefs, $ord_genomes_aref, $ord_kogs_aref ) = @_;
    $logger->info("Printing out the kog count table");
    open my $ALT, ">", "$out_dir/$group_name.count_tbl_rheat.txt";
    open my $OUT, ">", "$out_dir/$group_name.count_tbl.txt";
    
    #print the header
    print $ALT $group_name, "ID\tgenomeID\tcount\n";
    print $OUT $group_name, "ID\t", join("\t", @$ord_genomes_aref), "\n";
    
    #go through and print the information in the correct order
    foreach my $kog ( @$ord_kogs_aref ) {
        print $OUT "$kog";
        foreach my $geno ( @$ord_genomes_aref ) {            
            my $geno_href = $count_href_hrefs->{$geno};
            if ( !$geno_href->{$kog} ) {
                print $ALT "$kog\t$geno\t0\n";
                print $OUT "\t0";
                next;
            }
            print $OUT "\t", $geno_href->{$kog};
            print $ALT "$kog\t$geno\t", $geno_href->{$kog}, "\n";
        }
        print $OUT "\n";
    }
    
    close($ALT);
    close($OUT);
    return 1;
}


