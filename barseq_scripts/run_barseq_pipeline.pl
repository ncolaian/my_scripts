#! /usr/bin/env perl

#Need to have blat and r loaded

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

#My Variables

my $help = 0;
my $man  = 0;
my $organism_dir;
my $organism_name;
my $program_path;
my $set_dir_names;
my $bar_seq_data_path;
my $gc_file;
my $tsv;
my $rev;


# Read in variables from the command line
GetOptions ('man'  => \$man,
	    'help' => \$help,
	    'org_name=s' => \$organism_name,
	    'org_dir=s'  => \$organism_dir,
	    'prg_path=s' => \$program_path,
	    'bs_names=s' => \$set_dir_names,
	    'bs_path=s'  => \$bar_seq_data_path,
	    'tsv_name=s' => \$tsv,
	    'gc'         => \$gc_file,
			'rev'				 => \$rev,			
	   ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#Initiate Logger
my $logger = get_logger();

# MAIN #
#Load necessary modules
#system("module load r");
#system("module load blat");

#create directory and file objects 
my $main_do = dir($organism_dir);
my $set_name_fo = file($set_dir_names);

#run everything
#Will handle if the gbk file is passed, and change how the dir is set-up
if ( $gc_file ) {
  setup_w_gc( $main_do, $program_path);
}
else {
  setup_dir($main_do,$program_path);
}

#stall("set", $organism_dir);

local_bar_seq($organism_name, $organism_dir, $set_name_fo,
	      $bar_seq_data_path, $program_path);

stall("local", $organism_dir);

#r_bar_seq($organism_name, $organism_dir, $program_path);

### Subroutines that will occur in the main program ###

#Make sure that the directory contains a pool file, gff file, and a fna file
sub setup_dir {
  my ($dir_o,$prg_path) = @_;
  $logger->info("Setting up the directory for barseq programs");
  
  my $file_tracker = 0; #will keep track of what files are in the dir
  my $path = $dir_o->stringify();
  my $pool_file    = $path; #starting paths for files
  my $fna_file     = $path;
  my $gff_file     = $path;
  
  #1 = fna 2 = gff 3 = fna&gff 4 = pool 5 = pool&fna 6 = gff&pool 7 = all
  while (my $file = $dir_o->next) {
    next unless -f $file;
    next if -z $file;
    chomp $file;
    if ($file =~ m/\.gff/i) {
      $file_tracker += 2;
      $gff_file = $file;
    }
    elsif ($file =~ m/\.(fna|fasta)/i) {
      $file_tracker += 1;
      $fna_file = $file;
    }
    elsif ($file =~ m/\.pool/i) {
      $file_tracker += 4;
      $pool_file = $file;
    }
  }#While Loop

  #Gives appropriate error messages depending on what is missing
  if ($file_tracker != 7) {
    if ($file_tracker > 7) {
      croak("Too many files. Make sure there is only 1 fna, 1 gff, and 1 pool file\n");
    }
    elsif ($file_tracker == 6) {
      croak("Missing fna file\n");
    }
    elsif ($file_tracker == 5) {
      croak("Missing gff file\n");
    }
    elsif ($file_tracker == 4) {
      croak("Missing fna and gff file\n");
    }
    elsif ($file_tracker == 3) {
      croak("Missing pool file\n");
    }
    elsif ($file_tracker == 2) {
      croak("Missing fna and pool file\n");
    }
    elsif ($file_tracker == 1) {
      croak("Missing gff and pool file\n");
    }
    elsif ($file_tracker == 0) {
      croak("Check directory. None of the needed files are in the directory\n");
    }
  }# does not equal 7 if

  #check for tsv file
  if (!-e "$path/FEBA_BarSeq.tsv") {
    croak("Need to create a tsv file containing the barseq experiment metadata");
  }
  
  #Configure directory
  system("ln -s $pool_file $path/pool"); #creates a hard link
  system("bsub -o $path/lsf.out -J set perl $prg_path/SetupOrg.pl -gff $gff_file -fna $fna_file -out $path");
  mkdir("$path/html");
}

sub local_bar_seq {
  my ($org_name, $org_dir, $file_obj, $bs_path, $prg_path) = @_;
  $logger->info("Running local barseq with RunBarSeqLocal.pl");
  
  my $cmd = "bsub -M 20 -J local -e $org_dir/err.out -o $org_dir/lsf.out \"";

  my @seq_names = $file_obj->slurp(chomp=>1); #get bar seq run directory names
  foreach my $run_name (@seq_names) {
    my $path = $bs_path . "/" . $run_name;
    $cmd .= "perl $prg_path/RunBarSeqLocal.pl $org_dir $run_name $path";
		$cmd .= " -rev" if defined $rev;
		$cmd .= ";";
  }
  chop $cmd;
  $cmd .= "\"";
  
  $logger->debug($cmd);
  system($cmd);# will run all the local bar seqs in one run. 
}

sub r_bar_seq {
  my ($org_name, $org_dir, $prg_path) = @_;
  $logger->info("Running the r code - BarSeqR.pl");
  
  my $cmd = "bsub -o $org_dir/lsf.out perl $prg_path/BarSeqR.pl -org $org_name -indir $org_dir -outdir $org_dir/html";
  
  $logger->debug($cmd);
  system($cmd);
}

sub setup_w_gc {
  my ($dir_o,$prg_path) = @_;
  $logger->info("Setting up the directory for barseq programs");
  
  my $file_tracker = 0; #will keep track of what files are in the dir
  my $pool_file; 
  my $gc_f;
  my $path = $dir_o->stringify();
  
  #Go through the directory and find the files needed for setup
  while (my $file = $dir_o->next) {
    next unless -f $file;
    next if -z $file;
    chomp $file;
    
    if ($file =~ m/\.pool/i) {
      $file_tracker += 2;
      $pool_file = $file;
    }
    elsif ($file =~ m/\.GC/) {
      $file_tracker += 1;
      $gc_f = $file;
    }
  }
  
  #Produce error messages for files missing
  if ( $file_tracker != 3 ) {
    if ( $file_tracker == 0 ) {
      croak("$path does not include a genbank file (.gbk) or a pool file");
    }
    elsif ( $file_tracker == 1 ) {
      croak("$path does not include a pool file");
    }
    elsif ( $file_tracker == 2 ) {
      croak("$path does not include a .GC file");
    }
  }
  
  #check for tsv file
  if (!-e "$path/$tsv") {
    croak("Need to create a tsv file containing the barseq experiment metadata");
  }
  
  #Configure directory
  system("ln -s $pool_file $path/pool"); #creates a hard link
  mkdir("$path/html");
}

sub stall {
    my ( $job_name, $out_dir ) = @_;
    
    #Initiate the file
    my $cmd = "bjobs -J " . $job_name . " > $out_dir/ACTIVE_JOBS";
    system( $cmd );
    
    #Stall until file is empty
    while ( -s "$out_dir/ACTIVE_JOBS" ) {
        sleep 30;
        my $cmd = "bjobs -J " . $job_name . " > $out_dir/ACTIVE_JOBS";
        system( $cmd );
    }
    
    # Remove ACTIVE_JOBS file once the stalling is done
    `rm "$out_dir/ACTIVE_JOBS"`;
    
    return 1; 
}


__END__

=head1 Barseq Pipeline

This is a wrapper that will use a directory containing a pool file from a TnSeq experiment, a gff file, a fasta file, and a tsv file containing the BarSeq metadata to perform a barseq experiment. The output will be placed in a html directory within the organism directory.

The programs used are from FEBA created at the University of California. The source code is located at https://bitbucket.org/berkeleylab/feba

The programs being used are:
  SetupOrg.pl
  RunBarSeqLocal.pl
  BarSeqR.pl
  
Modules Needed To Be Loaded:
  R
  Blat

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

  NA

=head1 SYNOPSIS

  perl run_barseq_pipeline.pl -org_name -org_dir -prg_path -bs_names -bs_path [-gc]

  org_name => this is the name of the organism

  org_dir  => This is the path to the directory holding all the necessary files

  prg_path => path to the barseq bin (/proj/cdjones_lab/ncolaian/apps/feba/bin)

  bs_names => file with the names of the experiment data directories

  bs_path  => the path to the directory holding the barseq data

  tsv_name => The name of the tsv file

  gc      => This is a flag for if you have a gc file instead of an fna and gff file in the org dir

=head1 PARAMETERS

  org_name  =>  Organism Name     ->  String
  org_dir   =>  Directory Path    ->  String
  prg_path  =>  Path to FEBA bin  ->  String
  bs_names  =>  File Path         ->  String
  bs_path   =>  Directory Path    ->  String
  gbk       =>  Flag              ->  String  [Optional]
  
=head1  CONFIGURATION AND ENVIRONMENT
  
  Need to load Blat and R

=head1 Author

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=cut
