#! /usr/bin/env perl

# This object will be a scoring table that is created by passing two arrays
# containing a DNA sequence, the match scores

package Master_aln;

use strict;
use warnings;

use Getopt::Long;
use Carp;
use Readonly;
use List::Util qw(max);

#This new will create an object that can create a scoring table object
sub new{
  my ($class,%args) = @_;
 # my $scoring_table = Helpfxns::create_scoring_matrix($args{dnaseqi_ref},
				#		      $args{dnaseqj_ref});
  my $self = bless {
		   _scoring_table  => $args{scoring_table},
		   _dnaseqi_ref    => $args{dnaseqi_ref}
		                   || croak("Error: need 2 DNA sequences"),
		   _dnaseqj_ref    => $args{dnaseqj_ref}                                                           || croak("Error: need 2 DNA sequences"),
		   _match_score    => $args{match_score}
		                   || croak("Error: No match score"),
		   _mismatch_score => $args{mismatch_score}
		                   || croak("Error: No Mismatch score"),
		   _gap_score      => $args{gap_score}
		                   || croak("Error: No Gap Score"),
		  }, $class;

  $self -> set_scoring_table;

  return $self;
}
###########
# Getters #
###########
sub get_match_score    { $_[0] -> {_match_score}    }
sub get_mismatch_score { $_[0] -> {_mismatch_score} }
sub get_gap_score      { $_[0] -> {_gap_score}      }
sub get_dnaseqi_ref    { $_[0] -> {_dnaseqi_ref}    }
sub get_dnaseqj_ref    { $_[0] -> {_dnaseqj_ref}    }
sub get_scoring_matrix { $_[0] -> {_scoring_table}  }

sub get_scores         { return ($_[0] -> {_match_score   },
			         $_[0] -> {_mismatch_score},
			         $_[0] -> {_gap_score     })
		       }

sub get_imatrix_length { scalar(@{ ( $_[0] -> {_dnaseqi_ref} ) }) }
sub get_jmatrix_length { scalar(@{ ( $_[0] -> {_dnaseqj_ref} ) }) }

#This will return the i and j coordinate that represents the highest scoring
#array in the scoring matrix. It will also return the max scoring value
sub get_max_coordinates {
  my $self = $_[0];
  my $max_value = 0;
  my $i_coord = 0;
  my $j_coord = 0;
  
  my $imatrix_length = $self -> get_imatrix_length;
  my $jmatrix_length = $self -> get_jmatrix_length;
  my $scoring_matrix = $self -> {_scoring_table};

  for (my $i = 0; $i <= $imatrix_length; $i++){
    for (my $j = 0; $j <= $jmatrix_length; $j++){
      my $mayb_max_value = $scoring_matrix -> [$i][$j];
      if ($mayb_max_value > $max_value){
	$i_coord = $i;
	$j_coord = $j;
	$max_value = $mayb_max_value;
      }
      elsif ($mayb_max_value == $max_value){
	$i_coord = $i;
	$j_coord = $j;
      }
    }
  }
  return ($i_coord, $j_coord, $max_value);
}

###########
# Setters #
###########

#Cretes and sets the table in the object
sub set_scoring_table{
  my ($self) = @_;

  my $imatrix_length = $self -> get_imatrix_length;
  my $jmatrix_length = $self -> get_jmatrix_length;
  my $scoring_matrix = [];

  #create the matrix
  for (my $i = 0; $i <= $imatrix_length; $i++){
    for (my $j = 0; $j <= $jmatrix_length; $j++){
      #handle the first (0,0) array and the left and top sides (i=0 or j=0)
      if ($i == 0 && $j == 0){
	$scoring_matrix -> [$i][$j] = 0;
      }
      elsif ($i == 0){
	$scoring_matrix -> [$i][$j] = $j * $self -> {_gap_score};
      }
      elsif ($j == 0){
	$scoring_matrix -> [$i][$j] = $i * $self -> {_gap_score};
      }

      else{
	# Handle the rest of the scoring matrix by looking for the max value
	# taking into consideration the tile above, top-left, and left.
	
	my $diagonal = 0;
	my $delete = ($scoring_matrix -> [$i-1][$j]) + ($self -> {_gap_score});
	my $insert = ($scoring_matrix -> [$i][$j-1]) + ($self -> {_gap_score});
	# Determines if the diagonal is a match or mismatch
	if (($self -> {_dnaseqi_ref}) -> [$i-1] eq
	    ($self -> {_dnaseqj_ref}) -> [$j-1]){
	  $diagonal = ($scoring_matrix -> [$i-1][$j-1])
	  + ($self -> {_match_score});
	}
	else{
	  $diagonal = ($scoring_matrix -> [$i-1][$j-1])
	  + ($self -> {_mismatch_score});
	}
	$scoring_matrix -> [$i][$j] = max($diagonal,$delete,$insert);
      }
    }
  }
  $self -> {_scoring_table} = $scoring_matrix;
}

1;
