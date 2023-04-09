#!/usr/bin/env perl
# Copyright 2015 Oregon State University.
# All Rights Reserved. 
#
# This software program and documentation are copyrighted by OREGON STATE
# UNIVERSITY. The software program and documentation are supplied "as is", 
# without any accompanying services from the University. The University does not
# warrant that the operation of the program will be uninterrupted or error free. 
# The end-user understands that the program was developed for research 
# purposes and is advised not to rely exclusively on the program for any reason. 
#
# IN NO EVENT SHALL OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY 
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
# DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS 
# SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE OREGON STATE  
# UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY 
# STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED 
# HEREUNDER IS ON AN "AS IS" BASIS, AND OREGON STATE UNIVERSITY HAS 
# NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
# ENHANCEMENTS, OR MODIFICATIONS. 

##########################################################################
# Program : filterfastq.pl 
# Author  : Matthew Peterson
# Email   : matthew@cgrb.oregonstate.edu
# Created : 20150122
# Modified: 20150225
# Purpose : Filter out sequences that match a reference genome, e.g., PhiX
# from Illumina HiSeq/MiSeq gzipped FASTQ files (single or paired-end).
# Ideas for script design provided by Brent Kronmiller and Shawn O'Neil
##########################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Copy;
use Cwd 'abs_path';
# Unbuffer stdout for hash marks
$|=1;

# Binaries on infrastructure; MODIFY BINARY PATHS AS NEEDED
# NCBI BLAST
my $makeblastdb='/local/cluster/ncbi-blast-2.2.29+/bin/makeblastdb';
my $blastn='/local/cluster/ncbi-blast-2.2.29+/bin/blastn';
# FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
my $fastq_to_fasta='/local/cluster/bin/fastq_to_fasta';
# GNU split
my $split='/usr/bin/split';

# Test if binaries are executeable
if (! -x "$makeblastdb") {
    print "\nERROR: NCBI BLAST makeblastdb at '$makeblastdb' is not executeable.\n";
    print "       Modify \$makeblastdb variable in $0 to specify location.\n\n";
    exit(1);
}

if (! -x "$blastn") {
    print "\nERROR: NCBI BLAST blastn at '$blastn' is not executeable.\n";
    print "       Modify \$blastn variable in $0 to specify location.\n\n";
    exit(1);
}

if (! -x "$fastq_to_fasta") {
    print "\nERROR: FASTX-Toolkit fastq_to_fasta at '$fastq_to_fasta' is not executeable.\n";
    print "       Modify \$fastq_to_fasta variable in $0 to specify location.\n\n";
    exit(1);
}

if (! -x "$split") {
    print "\nERROR: GNU split at '$split' is not executeable.\n";
    print "       Modify \$split variable in $0 to specify location.\n\n";
    exit(1);
}
# Arguments to be passed
my $reference_genome; # FASTA reference genome to be blasted against
my $fastq_input_dir;  # Single directory containing gzipped FASTQ file(s) to be converted
my $fastq_output_dir; # Single directory to store filtered FASTQ and intermediary files
my $threads=1;        # Optional number of threads (default is 1)
my $evalue='1e15';    # Optional blastn evalue (default is: 1e-15)
my $quiet;            # Optional to not write script progress to standard output
my $help;             # Print out usage()

GetOptions(
    'reference|r=s' => \$reference_genome,
    'input|i=s' => \$fastq_input_dir,
    'output|o=s' => \$fastq_output_dir,
    'threads|t=i' => \$threads,
    'evalue|e=s' => \$evalue,
    'quiet|q' => \$quiet,
    'help|h' => \$help,
);

usage() if (defined $help);

# Pre-flight
usage() if ( not defined $reference_genome or not defined $fastq_input_dir or not defined $fastq_output_dir);

if (! -r "$reference_genome") {
    print "\nERROR: FASTA reference genome file '$reference_genome' is not readable.\n";
    usage();
}

if (! -d "$fastq_input_dir") {
    print "\nERROR: Directory '$fastq_input_dir' does not exist.\n";
    usage();
}

if (-d "$fastq_output_dir") {
    print "\nERROR: Directory '$fastq_output_dir' already exists; will not overwrite.\n";
    usage();
}

if (! -w (dirname "$fastq_output_dir")) {
    print "\nERROR: Cannot create output directory: '$fastq_output_dir'\n";
    usage();
}

# Check for number of gzipped FASTQ files
my @files = glob("$fastq_input_dir/*fastq.gz");
my $total_fastq = scalar(@files); # Count of entries in array
if ($total_fastq == 0) {
    print "\nERROR: No gzipped FASTQ files (*.fastq.gz) found in '$fastq_input_dir'\n";
    usage();
}

sub usage {
    print "\nUsage: $0 -r reference_genome.fa -i fastq_input_dir -o fastq_output_dir [-t threads] [-e evalue]\n\n";
    print "-r, --reference\n";
    print "       FASTA reference genome for sequences to be filtered against\n";
    print "-i, --input\n";
    print "       Directory containing the gzipped FASTQ file(s) to be filtered\n";
    print "       Sets of paired-end reads need to be in the same directory and are\n";
    print "       identified by their filenames containing either *_R1_* and *_R2_*\n";
    print "-o, --output\n";
    print "       Directory to store the filtered FASTQ and intermediary files\n";
    print "       Script will create output directory but not parents, e.g., mkdir -p\n";
    print "-t, --threads\n";
    print "       Optional number of threads for blastn to use (default is 1)\n";
    print "-e, --evalue\n";
    print "       Optional evalue for blastn to use (default is 1e-15)\n";
    print "-q, --quiet\n";
    print "       Do not write script progress to standard output\n";
    print "-h, --help\n";
    print "       This usage information.\n\n";
    exit(1);
}   

# Create output directory to hold filtered FASTQ files and intermediary files, i.e.,
# Copy of the reference genome, makeblastdb intermediary files, and blastn outputs
if (! defined $quiet) {
    print "Starting FASTQ filtering of $total_fastq files\n";
    print "Making output directory $fastq_output_dir\n";
}

unless(mkdir "$fastq_output_dir") {
    die "ERROR: Unable to create output directory: '$fastq_output_dir'";
}

my $fastq_output_dir_basename = basename $fastq_output_dir;
if (! defined $quiet) {
    print "Making intermediary directory $fastq_output_dir_basename/filterfastq\n";
}
unless(mkdir "$fastq_output_dir/filterfastq") {
    die "ERROR: Unable to create intermediary directory: '$fastq_output_dir/filterfastq'";
}

my $reference_genome_basename = basename $reference_genome;
# Copy reference genome to output directory
if (! defined $quiet) {
    print "Copying reference genome $reference_genome_basename to intermediary directory\n";
}
copy($reference_genome, "$fastq_output_dir/filterfastq/$reference_genome_basename") or die "ERROR: Reference genome '$reference_genome' cannot be copied.";

# Create BLAST database out of reference genome
if (! defined $quiet) {
    print "Creating BLAST database using $reference_genome_basename\n";
}
system("$makeblastdb -dbtype nucl -in $fastq_output_dir/filterfastq/$reference_genome_basename > /dev/null 2>&1");
my $makeblastdb_exit_code = $? >> 8;
if ( $makeblastdb_exit_code != 0)
{
    print "ERROR: $makeblastdb exited with code $makeblastdb_exit_code\n";
    exit(1);
}

# Itterate through all gzipped FASTQ files independent if there's a R2
# @files defined above during pre-flight checks
foreach my $file (@files) {
    my $file_basename = basename $file;
    # Specify -Q 33 to handle Illumina data using Sanger scores.
    # Specify -n to "keep sequences with unknown (N) nucleotides."
    # Header remains almost identical except for the first character where @ is transposed to >
    if (! defined $quiet) {
        print "Converting $file_basename from FASTQ to FASTA\n";
    }
    # sed statement will match any sequence starting with an N and repeating Ns to the 
    # end of the line and then remove the N sequence line and the header above it
    # Removing several sequences in a row contianing all Ns is necessary for blastn which can
    # potentially fail with "BLAST engine error" and an exit code of '3' due to low diversity
    # Use GNU split to create FASTA files containing 4 million reads each to ensure that 
    # BLAST finishes within the available memory on the system.
    # 
    # Since you change directory for GNU split before hand you need to determine the actual
    # path of the fastq_input_dir since supplying ./ as an input will no longer resolve correctly
    my $fastq_input_dir_abs_path = abs_path($fastq_input_dir);
    # Split on 8 million; each header with sequence is equal to two lines
    system("cd $fastq_output_dir/filterfastq/ ; gunzip -c $fastq_input_dir_abs_path/$file_basename | $fastq_to_fasta -n -Q33 | sed -n '/^N*\$/{s/.*//;x;d;};x;p;\${x;p;}' | sed '/^\$/d' | $split -l 8000000 - FASTASPLIT");
    my $fastq_to_fasta_exit_code = $? >> 8;
    if ( $fastq_to_fasta_exit_code != 0)
    {
        print "ERROR: $fastq_to_fasta filtering step exited with code $fastq_to_fasta_exit_code\n";
        exit(1);
    }

    # BLAST FASTA file against reference genome
    # stderr is surpressed due to blastn complaining about sequencings composed of Ns, e.g.,
    # "Warning: lcl|Query_9948 M01498:124:000000000-AB6JN:1:1101:9390:1655 1:N:0:1:
    #  Warning: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query
    #  sequence or its translation. Please verify the query sequence(s) and/or filtering options"
    # Itterate BLAST over each "FASTASPLIT" file containing a maximum of 4 million reads.
    if (! defined $quiet) {
        print "BLASTing $file_basename.fasta against $reference_genome_basename using $threads threads\n";
        # Check for number of intermediary FASTA files containing a maximum of 4 million reads each
        my @split_files = glob("$fastq_output_dir/filterfastq/FASTASPLIT*");
        my $total_splits = scalar(@split_files); # Count of entries in array
        print "$file_basename.fasta split into $total_splits intermediary FASTA files for processing\n";
    }
    
    # Create empty file to append to containing BLAST matches
    system("touch $fastq_output_dir/filterfastq/$file_basename.fasta.blastn.txt");
    my $touch_exit_code = $? >> 8;
    if ( $touch_exit_code != 0)
    {
        print "ERROR: touch attempting to create an empty file '$fastq_output_dir/filterfastq/$file_basename.fasta.blastn.txt' exited with code $touch_exit_code\n";
        exit(1);
    }

    # Itterate over each split file BLASTing it against the reference and append the blastn.txt output file
    my @split_files = glob("$fastq_output_dir/filterfastq/FASTASPLIT*");
    foreach my $split_file (@split_files) {
        my $split_file_basename = basename $split_file;
        # Print out a hash mark to indicate progress at the start of each blastn job
        if (! defined $quiet) {
            print ".";
        }
        # Start BLAST on intermediary file
        system("$blastn -db $fastq_output_dir/filterfastq/$reference_genome_basename -query $fastq_output_dir/filterfastq/$split_file_basename -max_target_seqs 1 -evalue $evalue -outfmt 6 -num_threads $threads 2> /dev/null >> $fastq_output_dir/filterfastq/$file_basename.fasta.blastn.txt");
        my $blastn_exit_code = $? >> 8;
        if ( $blastn_exit_code != 0)
        {
            print "ERROR: $blastn exited with code $blastn_exit_code\n";
            if($blastn_exit_code == 3) {
                print "       Exit code 3 may indicate numerous 'low-complexity' sequences";
                print "       e.g., sequences with too many Ns. See also:";
                print "       http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#LCR";
            }
            if($blastn_exit_code == 137) {
                print "       Exit code 137 indicates that $blastn ran out of memory!";
                print "       Consider reducing the number of threads or running with more RAM.";
            }
            exit(1);
        }
        # Remove intermediary SPLIT file as it's no longer needed
        unlink "$fastq_output_dir/filterfastq/$split_file_basename"
    }
    # Finished BLASTing against SPLIT files, print newline
    if (! defined $quiet) {
        print "\nFinished BLASTing $file_basename.fasta against $reference_genome_basename\n";
    }
    
} 

# @files defined above during pre-flight checks
# For each FASTQ file to be filtered
foreach my $file (@files) {
    my $file_basename = basename $file;
    if (! defined $quiet) {
        print "Reading original $file\n";
    }

    # Create hash based upon BLAST matches for both single and paired-end reads
    my %blast_hash = ();
    # If there are pair-end reads you need to load both sets of BLAST matches
    my $blast_file_wildcard = $file_basename;
    # Pair-end reads are denoted by containing either _R1_ or _R2_ in the filename
    $blast_file_wildcard =~ s/_R[12]_/_R{1,2}_/;
    my @blast_files = glob("$fastq_output_dir/filterfastq/$blast_file_wildcard.fasta.blastn.txt");

    # Load each BLAST file into hash
    foreach my $blast_file (@blast_files) {
        my $blast_file_basename = basename $blast_file;
        # There may not be an R2
        if (-r "$blast_file") {
            if (! defined $quiet) {
                print "Reading BLAST hits from $blast_file_basename\n";
            }
            open FH, "$blast_file" or die "ERROR: Unable to open $blast_file\n";
            while (my $line = <FH>) {
                # Split on BLAST output to grab first column and assign to $fasta_header_id
                my ($fasta_header_id) = split /\t/, $line;
                # Prefix with FASTQ '@' header; key no value
                $blast_hash{'@' .  $fasta_header_id} = undef;
            }
        }
    }

    my $total_blast_hits = keys %blast_hash;
    if (! defined $quiet) {
        print "Found a total of $total_blast_hits BLAST hits\n"; 
    }

    # Open gzipped output file for writing
    my $output_fastq = "$fastq_output_dir/$file_basename";
    if (! defined $quiet) {
        print "Writing filtered output to $output_fastq\n";
    }
    open OFH, "| /bin/gzip -c > $output_fastq" or die "ERROR: Unable to open $output_fastq\n";
    # Itterate over gzip file outputing reads that did not match the BLAST hits
    open IFH, "/bin/gunzip -c $file | " or die "ERROR: Unable ot open $file\n";
    my $fastq_count = 0;
    my $fastq_header_count = 0;
    while(my $line  = <IFH>) {
        # Count number of sequences
        $fastq_count++;
        # Header
        my $fastq_header = $line;
        # Sequence
        $line = <IFH>;
        my $fastq_sequence = $line;
        # Third row '+'
        $line = <IFH>;
        my $fastq_plus = $line;
        # Quality
        $line = <IFH>;
        my $fastq_quality = $line; 

        # Split on fasta_header to grab first word and assign to $fastq_header_id
        my ($fastq_header_id) = split / /, $fastq_header;
        # If the header is not found in the hash then print the sequence
        if (! exists $blast_hash{$fastq_header_id}) {
            print OFH $fastq_header;
            $fastq_header_count++;
            # Print hash marks every 1 million reads
            if (($fastq_header_count%1000000)==0) {
                if (! defined $quiet) {
                    print '.';
                }
            }
            print OFH $fastq_sequence;
            print OFH $fastq_plus;
            print OFH $fastq_quality;
        }
    } 
    close IFH;
    close OFH;
    my $filtered_percentage = 0;
    if ($fastq_count != 0) {
        $filtered_percentage = sprintf("%.3f", ((($fastq_count - $fastq_header_count)/$fastq_count) * 100));
    }
    if (! defined $quiet) {
        print "\nA total of $fastq_header_count out of $fastq_count sequences written ($filtered_percentage% filtered out)\n";
    }
}

# Copyright 2015 Oregon State University.
# All Rights Reserved. 
# 
# petersm3@cgrb.oregonstate.edu
