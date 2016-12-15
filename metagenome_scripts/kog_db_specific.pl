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
my $fungal_dir;
my $genome_names;
my $full_kog_file;
my $out;

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'full_kog=s' => \$full_kog_file,
            'fung_dir=s' => \$fungal_dir,
            'nm_file=s'  => \$genome_names,
            'out=s'      => \$out,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Subroutines that will occur in the main program

sub add_kogs_to_full_kog_list {
    my ($id, $dir, $kog_list_href) = @_;
    my $file_path = "$dir/$id/$id.kog.tab";
    #create file obj
    my $file_obj = file($file_path);
    my @kog_line_arefs = $file_obj->slurp( chomp=>1, split=>"\t" );
    
    #find id, class, and group cols
    my $col_nums = find_id_class_group_cols( $kog_line_arefs[0] );
    #add kogs to href
    foreach my $kog_line_aref ( @kog_line_arefs ) {
        if ( $kog_line_aref->[0] =~ qr/transcript/i ) {
            next;
        }
        if ( !exists $kog_list_href->{$kog_line_aref->[$col_nums->[0]]} ) {
            $kog_list_href->{$kog_line_aref->[$col_nums->[0]]} =
            [ $kog_line_aref->[$col_nums->[1]], $kog_line_aref->[$col_nums->[2]] ];
        }
        else {
            next;
        }
    }
    
    return 1;
}

sub find_id_class_group_cols {
    my ($line_aref) = @_;
    my $col_nums = [0,0,0];
    
    for (my $i = 0; $i < scalar( @$line_aref ); $i++) {
        if ( $line_aref->[$i] eq "kogid" ) {
            $col_nums->[0] = $i;
        }
        elsif ( $line_aref->[$i] eq "kogClass" ) {
            $col_nums->[1] = $i;
        }
        elsif ( $line_aref->[$i] eq "kogGroup" ) {
            $col_nums->[2] = $i;
        }
    }
    
    return ( $col_nums );
}

sub print_directory_specific_kog_list {
    my ( $kog_href, $kog_aref, $outfile ) = @_;
    open( my $fh, ">", $outfile );
    #print header
    print $fh "kogID\tkogName\tkogClass\tkogGroup\n";
    my @kogs_removed;
    
    foreach my $ordered_kog_aref ( @$kog_aref ) {
        my $id = $ordered_kog_aref->[0];
        if ( exists $kog_href->{ $id } ) {
            print $fh join( "\t", @$ordered_kog_aref), "\t";
            print $fh join ( "\t", @{$kog_href->{ $id }} ), "\n";
        }
        else {
            push( @kogs_removed, $id );
        }
    }
    #Tell the user what kogs were removed
    if ( scalar(@kogs_removed) > 0 ) {
        my $statement = join( ", ", @kogs_removed );
        carp "Kogs Removed: $statement";
    }
    
    close($fh);
    return 1;
}

# MAIN #
my $full_kog_list = {}; #holds kogids as keys and then an aref w/ kog class and group
my @kog_order_list; # Contains KogId and Name

#Set up File objects
my $full_kog_fo = file($full_kog_file);
my $genome_names_fo = file($genome_names);

#create the kog order array
my @full_kog_array_aref = $full_kog_fo->slurp( chomp=>1, split=>qr/\t/ );

foreach my $line_aref ( @full_kog_array_aref ) {
    if ($line_aref->[1] eq "KOG Name") {
        next;
    }
    push @kog_order_list, $line_aref;
}

# Read in file names and open and read kog files
my @names = $genome_names_fo->slurp( chomp=>1 );

foreach my $genome_id ( @names ) {
    add_kogs_to_full_kog_list( $genome_id, $fungal_dir, $full_kog_list );
}

print_directory_specific_kog_list( $full_kog_list, \@kog_order_list, $out );

__END__

=head1 

=head1 VERSION

This documentation refers to kog_db_specific.pl 0.0.1

=head1 INCLUDED MODULES



=head1 INHERIT

    NA
    
=head1 SYNOPSIS

    
    
=head1 PARAMETERS
    
     

=head1 CONFIGURATION AND ENVIRONMENT

    NA

=head1 DEPENDENCIES

    
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

    No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 Author

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2016, Nicholas Colaianni
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
