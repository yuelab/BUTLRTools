#! /usr/local/bin/perl -w

#Yanli Wang
#3dgenome.browser@gmail.com
#Version 1.3
#Converts the matrix output from HOMER to separate interchromosomal and intrachromosomal matrix files, as well as a list file
#Both would serve as input for BUTLR file conversion.

use strict;
use warnings;

#Find the path of the required modules
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/.");


use Getopt::Long qw(GetOptions);
use File::Basename;
use Butlr;

use constant COL_START_INDEX => 2;

my $version = "1.3";
my $homer_matrix;
my $genome_size_filename;
my $output_prefix;

#Print the usage when the required parameters are not provided
sub get_use
{
    return "BUTLR conversion tool version $version\n" . 
    "Converts HOMER one-matrix output to separate matrices and list file as inputs to create BUTLR files\n\n".
    "Usage: perl $0 <REQUIRED> <OPTIONAL>\n\n\t".
    "-m <name of homer matrix file> [REQUIRED]\n\n\t".
    "-g <genome size file/*.chrom.sizes> [REQUIRED]\n\n\t".
    "-o <output file prefix> [OPTIONAL]\n\n\n";
}

GetOptions(
    'matrix|m=s' => \$homer_matrix,
    'genome|g=s' => \$genome_size_filename,
    'output|o=s' => \$output_prefix,
) or die get_use();

unless (defined($homer_matrix) &&
        defined($genome_size_filename))
{
    die get_use();
}

if ( ! -e $genome_size_filename) {die "$! ($genome_size_filename)\n"};
my %genome_size;
open(FILE, $genome_size_filename) or die("Error opening ($genome_size_filename)\n");
    while (my $line = <FILE>)
    {
        chomp $line;
        my ($chr, $size) = split(/\s+/, $line);
        $genome_size{$chr} = $size;
    }
close(FILE);

if ( ! -e $homer_matrix) {die "$! ($homer_matrix)\n"};
open(FILE, $homer_matrix) or die("Error opening ($homer_matrix)\n");

my $prev_chrom1 = '';
my @col_title = ();
my $header_flag = 1;
my $new_row_flag = 0;

my ($basename, $dir, $ext) = fileparse($homer_matrix, qr/\.[^.]*/);

my $user_input_name_flag = 1;
if ( not defined $output_prefix )
{
    $user_input_name_flag = 0;   
}
else
{
    if ( -d $output_prefix )
    {
        $user_input_name_flag = 0;
    }
}

my $list_filename;
if ( !$user_input_name_flag )
{
    $list_filename = "$dir/matrix.list";
}
else
{
    $list_filename = "$output_prefix.list";
}
open LISTFILE, ">$list_filename" or die;

#With chrom1 as the name of each chromosome listed vertically in column 1 and chrom2 as the name of each chromosome
#listed horizontally in row 1. 
#Exact the matrix when the chromosome name changes, indicating the end of each interchromosomal/intrachromosomal 
#chromosomes
#To avoid redundancy, i.e. chr1 x chr2 and chr2 x chr1 store the same information, each row will always be the more
#From version 1.3, resulting matrices have # columns ( awk -F "\t" {print NF} [matrix] | uniq ) >= # rows ( wc -l [matrix] ) instead of the other way around.
while (my $line = <FILE>)
{
    chomp $line;
    my @recs = split(/\t/, $line);

    #One line of header, which contains the name of chromosome and start location of bin
    if ( $header_flag )
    {
        @col_title = @recs;
        $header_flag = 0;
        next;
    }

    my $chrom1 = substr( $recs[0], 0, rindex($recs[0], "-") );
    if ( !exists($genome_size{$chrom1}) )
    {
        next;   
    }

    if ( $prev_chrom1 ne $chrom1 )
    {
        $new_row_flag = 1;
        $prev_chrom1  = $chrom1;
    }

    my $new_col_flag = 0;
    my $prev_chrom2  = '';
    my $prev_index   = COL_START_INDEX;
    my $curr_index   = 0;

    #Traverse each column
    for my $i (COL_START_INDEX .. $#recs) 
    {
        my $chrom2 = substr( $col_title[$i], 0, rindex($col_title[$i], "-") );

        my $output_matrix_filename;
        #Note the two chromosome/scaffold names here and tentatively determine the output filename
        if ( !$user_input_name_flag )
        {
            $output_matrix_filename = "$dir$chrom1.$chrom2.matrix";
        }
        else
        {
            $output_matrix_filename = "$output_prefix.$chrom1.$chrom2.matrix";
        }

        if ( $prev_chrom2 ne $chrom2 )
        {
            $new_col_flag = 1;
            $curr_index   = $i;

            if ($prev_chrom2 ne '')
            {
                if ( !$user_input_name_flag )
                {
                    $output_matrix_filename = "$dir$chrom1.$prev_chrom2.matrix";
                }
                else
                {
                    $output_matrix_filename = "$output_prefix.$chrom1.$prev_chrom2.matrix";
                }
            }

            if ( $new_row_flag && exists($genome_size{$prev_chrom2}) && !is_chrom1_ahead(\%genome_size, $chrom1, $prev_chrom2) )
            {
                open MATXFILE, ">$output_matrix_filename" or die;
                close MATXFILE;
            }

            if ($prev_index < $curr_index && exists($genome_size{$prev_chrom2}) && !is_chrom1_ahead(\%genome_size, $chrom1, $prev_chrom2) )
            {
                open MATXFILE, ">>$output_matrix_filename" or die;
                print MATXFILE join( "\t", @recs[$prev_index .. $curr_index-1] ) . "\n";
                close MATXFILE;
            }
            $prev_index = $curr_index;
            $prev_chrom2  = $chrom2;
        }

        if ( $new_row_flag && $new_col_flag )
        {
            $new_col_flag = 0;

            if ( exists($genome_size{$chrom2}) && !is_chrom1_ahead(\%genome_size, $chrom1, $chrom2) )
            {
                if ( !$user_input_name_flag )
                {
                    $output_matrix_filename = "$dir$chrom1.$prev_chrom2.matrix";
                }
                else
                {
                    $output_matrix_filename = "$output_prefix.$chrom1.$prev_chrom2.matrix";
                }
                print LISTFILE $chrom1 . "\t" . $chrom2 . "\t" . "$output_matrix_filename" . "\n";
            }
        }
    }

    #Print the matrix when end of the matrix is denoted with column/file end, instead of chromosome name change.
    if ( exists($genome_size{$prev_chrom2}) && !is_chrom1_ahead(\%genome_size, $chrom1, $prev_chrom2) )
    {
        my $output_matrix_filename;
        if ( !$user_input_name_flag )
        {
            $output_matrix_filename = "$dir$chrom1.$prev_chrom2.matrix";
        }
        else
        {
            $output_matrix_filename = "$output_prefix.$chrom1.$prev_chrom2.matrix";
        }

        if ( $new_row_flag )
        {
            open MATXFILE, ">$output_matrix_filename" or die;
            close MATXFILE;
        }
        open MATXFILE, ">>$output_matrix_filename" or die;
        if ($prev_index == $#recs) #one column of info
        {
            print MATXFILE $recs[$#recs]  . "\n";
        }
        else #multiple columns of info
        {
            print MATXFILE join( "\t", @recs[$prev_index .. $#recs] ) . "\n";
        }
        close MATXFILE;
    }
    $new_row_flag = 0;
}

close LISTFILE;
close FILE;
