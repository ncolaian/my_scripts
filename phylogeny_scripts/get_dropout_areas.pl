#! usr/bin/evn perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Data::Dumper;
use Statistics::Descriptive;


# My Variables
my $help = 0;
my $man = 0;
my $name_file;
my $gene_name; # Will hold the param file passed into the driver
my $out_dir; #place to put the ordered ids and meta_file

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'name_file=s'       => \$name_file,
            'out_dir|o=s'       => \$out_dir,
            'gene_name|gn=s'    => \$gene_name,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Subroutines that will occur in the main program

sub find_mean_std {
    my ( $files_aref, $mean_href, $std_href, $names_aref) = @_;
    my @count_to_position_hrefs;
    
    for (my $i = 0; $i < scalar @$files_aref; $i++) {
        my @counts;
        my $count_to_pos_href = {};
        
        my $file_o = file( $files_aref->[$i] );
        my @lines = $file_o->slurp( chomp=>1, split=>qr/,/ );
        
        #get counts and push a holder of position and count
        foreach my $line_aref ( @lines ) {
            if ( exists $count_to_pos_href->{$line_aref->[1]} ) {
                push @{ $count_to_pos_href->{$line_aref->[1]} }, $line_aref->[0];
            }
            else {
                $count_to_pos_href->{$line_aref->[1]} = [$line_aref->[0]];
            }
            push @counts, $line_aref->[1];
        }
        #Calculate mean and std
        my $stat_obj = Statistics::Descriptive::Full->new();
        $stat_obj->add_data(@counts);
        my $mean = $stat_obj->mean();
        my $std = $stat_obj->standard_deviation();
        
        $mean_href->{$names_aref->[$i]} = $mean;
        $std_href->{$names_aref->[$i]} = $std;
        push @count_to_position_hrefs, $count_to_pos_href;
    }
    return \@count_to_position_hrefs;
}

sub print_indiv_file_stats {
    my ( $name, $mean, $std, $out, $count_href ) = @_;
    my $excluded_aref = $count_href->{0};
    
    if ( $mean - ($std*3) >= 1 ) {
        my $min = $mean - ($std*3);
        $excluded_aref = $count_href->{0};
        for (my $i = 1; $i < $min; $i++) {
            foreach my $pos ( @{ $count_href->{$i} } ) {
                push @$excluded_aref, $pos;
            }
        }
    }
    my $drops = scalar @$excluded_aref;
    my @sorted_excluded = sort { $a <=> $b } @$excluded_aref;
    open (my $fh, ">>", $out);
    print $fh "$name\tmean=$mean\tstd=$std\tdrops=$drops\n";
    print $fh join( ",", @sorted_excluded), "\n\n";
    close($fh);
}

# MAIN #
my %means;
my %std;
my @paths;
my $file_dir = "/proj/cdjones_lab/projects/waddell/coverage/results";
my $file_extension = ".asm_genes.pileup."; #add the gene name and .csv to end

#get names from file
my $name_fo = file( $name_file );
my @names = $name_fo->slurp( chomp=>1 );
foreach my $name ( @names ) {
    my $file_path = "$file_dir/$name/$name" . $file_extension . "$gene_name.csv";
    push @paths, $file_path;
}

my $count_to_pos_hrefs_aref = find_mean_std(\@paths, \%means, \%std, \@names);
#Own sanity check
print Dumper(%means) , "\n";
if ( !%means || !%std ) {
    print STDERR "The mean and std's were not saved";
}

#print stuff out
my $output = "$out_dir/$gene_name.dropout_pos.txt";
for ( my $i = 0; $i < scalar(@names); $i++) {
    print_indiv_file_stats( $names[$i], $means{$names[$i]}, $std{$names[$i]},
                            $output, $count_to_pos_hrefs_aref->[$i] );
}




