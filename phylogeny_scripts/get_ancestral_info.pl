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
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::NodeI;



# My Variables
my $help = 0;
my $man = 0;
my $newick_file;
my $output;

# Read in the variable from the command line
GetOptions ( 'man'  =>  \$man,
            'help'  =>  \$help,
            'tree|t=s'      => \$newick_file,
            'output|o=s'    => \$output,
            ) || die("There was an error in the command line arguements \n");

# Pod Usage for the manual and help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

## MAIN ##
$logger->info( "Loading Tree File" );
my $input = Bio::TreeIO->new( -file => $newick_file,
                              -format => "newick" );
$logger->info( "creating a tree object" );
my $tree = $input->next_tree;

$logger->info( "getting nodes" );
my $node = $tree->get_root_node;

my @nodes = $node->each_Descendent;
my $an_num = 1;

my $i = 1;
my $level = 1;
while ( $i == 1 ) {
    print "--$level:\n";
    my $non_leaf = 0;
    my @new_nodes;
    foreach my $n ( @nodes ) {
        if ( $n->id eq "" ) {
            $n->id("Ancestor$an_num");
            $an_num++;
        }
        
        if ( $n->is_Leaf ) {
            print $n->id, "\tLeaf\n";
        }
        else {
            $non_leaf++;
            my @id_array;
            my @mm_array = $n->each_Descendent;
            foreach my $part ( @mm_array ) {
                if ( $part->id eq "" ) {
                    $part->id("Ancestor$an_num");
                    $an_num++;
                }
                push @id_array,$part->id;
            }
            print $n->id, "\t", join ", ", @id_array, "\n";
            push @new_nodes, @mm_array;
        }
    }
    if ($non_leaf == 0) {
        $i++;
    }
    $level++;
    @nodes = @new_nodes;
}



