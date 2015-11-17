#! /usr/local/bin/perl -w

#Yanli Wang
#yanly.wang@gmail.com
#Version 1.2
#Code to extract Hi-C matrix (contact matrices in multiple files) from the Binary Upper TrianguLer MatRix (BUTLR) format

use strict;
use warnings;

#Find the path of the required modules
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/.");

#Required modules
use Getopt::Long qw(GetOptions);
use File::Basename;
use List::Util qw(any first);
use Butlr;

my $version = "1.2";

my $butlr_filename;
my $output_prefix;
my $output_filename;
my $location_str;
my $bin_str;

my $chrom1_name  = "";
my $chrom2_name  = "";


my $UNINIT_VALUE     = -0.123456789;
my $chrom1_start     = $UNINIT_VALUE;
my $chrom1_endin     = $UNINIT_VALUE;
my $chrom1_start_bin = $UNINIT_VALUE;
my $chrom1_endin_bin = $UNINIT_VALUE;

my $chrom2_start     = $UNINIT_VALUE;
my $chrom2_endin     = $UNINIT_VALUE;
my $chrom2_start_bin = $UNINIT_VALUE;
my $chrom2_endin_bin = $UNINIT_VALUE;

my $DATA_BYTE_SIZE = 4;
my $LOCN_BYTE_SIZE = 8;

sub get_use
{
    return  "BUTLR conversion tool version $version\n".
            "Usage: perl $0 <REQUIRED> <OPTIONAL>\n\t".
            "-i <input butlr filename> [REQUIRED]\n\t".
            "-l <location, 0-based, format: chromosome1[:start1-end1][,chromosome2[:start2-end2]]> [OPTIONAL]\n\t".
            "\tExtract the matrix based on the location. Examples: chr1:0-100000 (both row and col of matrix) or\n\t".
            "\tchr1:0-100000,249150621-249250621 (row and col of same chromosome/intrachromosomal) or chr1 (entire chromosome) \n\t".
            "\tor chr1:0-100000,chr2:0-100000 or chr1,chr2 (entire interchromosomal) interaction\n\t".
            "\t*order matters: first chromosome is the row while second is the column.\n\t".
            "-b <bin (= location / resolution), 0-based, use bin instead of absolute position> [OPTIONAL]\n\t".
            "-o <output prefix> [OPTIONAL]\n\n".
            "\n\n";
            #"\tIf neither location or bin is provided, all chromosomes/scaffold are written to the output file\n\n";
}

GetOptions(
    'input|i=s' => \$butlr_filename,
    'location|l=s' => \$location_str,
    'bin|b=s' => \$bin_str,
    'output|o=s' => \$output_prefix,
) or die get_use();

unless (defined($butlr_filename))
{
    die get_use();
}

if (defined($location_str) && defined($bin_str))
{
    print STDERR "Both location and bin are defined. Please select one or the other.\n";
    die get_use();
}

#Input:
#   1) string: "chrom[:start-end]"
#Ouput:
#   array: (chrom, start, end)
sub extract_location
{
    my @loc_temp = split(":", $_[0]);
    my $c = $loc_temp[0];

    if    ( scalar @loc_temp == 1 )
    {
        return ($c, $UNINIT_VALUE, $UNINIT_VALUE);
    }
    elsif ( scalar @loc_temp == 2 )
    {
        my @rec_temp  = split("-", $loc_temp[1]);
        if    (scalar @rec_temp == 1 )
        {
            return ($c, $rec_temp[0], $UNINIT_VALUE);
        }
        elsif (scalar @rec_temp == 2)
        {
            return ($c, $rec_temp[0], $rec_temp[1]);
        }
        else #Inappropriate parsing of ':'
        {
            print STDERR "Error parsing the chromosome location information.\n\n";
            die get_use();
        }
    }
    else
    {
        print STDERR "Error parsing the chromosome location information.\n\n";
        die get_use();   
    }
}

my ($basename, $dir, $ext) = fileparse($butlr_filename, qr/\.[^.]*/);
if (not defined $output_prefix)
{
    $output_prefix = "$dir$basename";
}

#Input:  
#   1) Pointer to BUTLR file
#   2) Number of bytes to read
#Subroutine:
#   Read binary file into buffer
#Output:
#   Binary String
sub read_bytes
{
    my $butlr_fp = $_[0];
    my $byte_num = $_[1];
    my $buffer;
    my $n;
    if ( $n = read( $butlr_fp, $buffer, $byte_num ) )
    {
        if ( $n != $byte_num )
        {
            die "[Error] File reading has yielded unexpected results (incorrect file size?).\n";
        }
    }
    else
    {
        die "[Error] File reading has failed.\n";   
    }
    return $buffer;  
}
#Input:  
#   1) Pointer to BUTLR file
#Subroutine:
#   Reads binary files as characters and convert to string
#Output:
#   String
sub read_chars
{
    my $butlr_fp =$_[0];
    my $c;
    my $string = '';
    while ( my $n = read $butlr_fp, $c, 1 )
    {
        my @clist;
        push @clist, unpack('C', $c);
        if ( unpack('C', $c) == 0 )
        {
            last;
        }
        $string .= convert_from_bin_to_char(\@clist);
    }
    return $string;
}
#Input:  
#   1) Pointer to BUTLR file
#   2) Location 
#   3) Row
#   4) Column Start
#   5) Column End
#Subroutine:
#   Read butlr file to extract matrix
#Output:
#   Array of values within the provided row and column of the interactions matrix
sub get_values
{
    my $butlr_fp =  $_[0];
    my $location =  $_[1];
    my $row_offset= $_[2]; #i
    my $col_start = $_[3]; #j
    my $col_endin = $_[4];

    if ($row_offset > $col_start && $col_start == $col_endin)
    {
        my $temp = $col_start;
        $col_start  = $row_offset;
        $col_endin  = $row_offset;
        $row_offset = $col_start;
    }

    seek( $butlr_fp, $location + $LOCN_BYTE_SIZE * $row_offset, 0 );
    my $row_locn_start = unpack('Q', read_bytes( $butlr_fp, $LOCN_BYTE_SIZE ));
    seek( $butlr_fp, $location + $LOCN_BYTE_SIZE * ($row_offset+1), 0 );
    my $row_locn_endin = unpack('Q', read_bytes( $butlr_fp, $LOCN_BYTE_SIZE ));

    my %index_to_value;
    for (my $i = $row_locn_start; $i < $row_locn_endin; $i+=$DATA_BYTE_SIZE*2)
    {
        seek( $butlr_fp, $i, 0 );
        my $l = unpack('L',  read_bytes($butlr_fp, $DATA_BYTE_SIZE));

        if ($l >= $col_start && $l <= $col_endin)
        {
            my $v = unpack('f<', read_bytes($butlr_fp, $DATA_BYTE_SIZE));
            $index_to_value{$l} = $v;
        }
        if ($l > $col_endin) { last; }
    }
    my @list;
    for (my $i = $col_start; $i <= $col_endin; $i++)
    {
        if (exists $index_to_value{$i}) { push @list, $index_to_value{$i}; }
        else { push @list, 0.0; }
    }
    return @list;
}

if ( ! -e $butlr_filename ) { die "$! ($butlr_filename)\n" };
open(my $butlr_fp, $butlr_filename) or die "Error opening ($butlr_filename)\n";

my $header_size = unpack('L', read_bytes($butlr_fp, $DATA_BYTE_SIZE));

my @clist = unpack('C*', read_bytes($butlr_fp, 16));
print STDERR " BUTLR Converter Version: " . convert_from_bin_to_char(\@clist) . "\n";

print STDERR " Size of header: " . $header_size . "\n";

my $intra_locn = unpack('L', read_bytes($butlr_fp, $DATA_BYTE_SIZE));
print STDERR " Location of chr information: "  . $intra_locn . "\n";

my $inter_locn = unpack('L', read_bytes($butlr_fp, $DATA_BYTE_SIZE));
print STDERR " Location of interchromosomal: " . $inter_locn . "\n";

my $assembly = read_chars( $butlr_fp );
print STDERR " Assembly: " . $assembly . "\n";

my $res = unpack('L', read_bytes($butlr_fp, $DATA_BYTE_SIZE));
print STDERR " Resolution: " . $res . "-bp\n";

seek( $butlr_fp, 20, 1 );

my @sorted_chr_list;
my %chr_to_size;
my %chr_to_locn;
print STDERR "chromosome/scaffold\tsize\tlocation\n";

my $curr_byte = $intra_locn;
my $intra_end = $inter_locn;
if (! $intra_end) { $intra_end = $header_size; }
while ($curr_byte < $intra_end)
{
    my $chr  = read_chars( $butlr_fp );
    push @sorted_chr_list, $chr;
    my $size = unpack('L', read_bytes($butlr_fp, $DATA_BYTE_SIZE));
    $chr_to_size{ $chr } = $size;
    my $location = unpack('Q', read_bytes($butlr_fp, $LOCN_BYTE_SIZE));
    $chr_to_locn{ $chr } = $location;
    $curr_byte += length($chr) + 1 + $DATA_BYTE_SIZE + $LOCN_BYTE_SIZE;
    print STDERR "\t$chr\t$size\t$location\n"
}


if ( $inter_locn )
{
    print STDERR "chromosome/scaffold1\tchromosome/scaffold2\tlocation\n";
while ( $curr_byte < $header_size )
{
    my $key  = read_chars( $butlr_fp );
    my $location = unpack('Q', read_bytes($butlr_fp, $LOCN_BYTE_SIZE));
    $curr_byte += length($key) + 1 + $LOCN_BYTE_SIZE;
    print STDERR "\t$key\t$location\n";
}
}

my $chrom1_rowcol_num;
my $chrom2_rowcol_num;

if (defined($location_str))
{   
    my @chr_temp = split(",", $location_str);
    ($chrom1_name, $chrom1_start, $chrom1_endin) = extract_location( $chr_temp[0] );
    
    if    ( scalar @chr_temp == 1 )
    {   
        ($chrom2_name, $chrom2_start, $chrom2_endin) = ($chrom1_name, $chrom1_start, $chrom1_endin);
    }
    elsif ( scalar @chr_temp == 2 )
    {   
        ($chrom2_name, $chrom2_start, $chrom2_endin) = extract_location( $chr_temp[1] );
    }
    else
    {   
        print STDERR "Error parsing the chromosome location information.\n\n";
        die get_use();
    }

    if ( !exists $chr_to_size{ $chrom1_name } )
    {
        die "Chromosome/Scaffold $chrom1_name not in the file.\n";
    }
    if ( !exists $chr_to_size{ $chrom2_name } )
    {
        die "Chromosome/Scaffold $chrom2_name not in the file.\n";   
    }
    $chrom1_rowcol_num = int( $chr_to_size{ $chrom1_name } / $res ) + 1;
    $chrom2_rowcol_num = int( $chr_to_size{ $chrom2_name } / $res ) + 1;

    if ( $chrom1_start == $UNINIT_VALUE ) { $chrom1_start = 0; }
    if ( $chrom2_start == $UNINIT_VALUE ) { $chrom2_start = 0; }
    if ( $chrom1_endin == $UNINIT_VALUE ) { $chrom1_endin = $chr_to_size{ $chrom1_name }; }
    if ( $chrom2_endin == $UNINIT_VALUE ) { $chrom2_endin = $chr_to_size{ $chrom2_name }; }

    $chrom1_start_bin = int($chrom1_start / $res);
    $chrom1_endin_bin = int($chrom1_endin / $res);
    $chrom2_start_bin = int($chrom2_start / $res);
    $chrom2_endin_bin = int($chrom2_endin / $res);
}

if (defined($bin_str))
{
    my @chr_temp = split(",", $bin_str);
    ($chrom1_name, $chrom1_start_bin, $chrom1_endin_bin) = extract_location( $chr_temp[0] );
    
    if    ( scalar @chr_temp == 1 )
    {   
        ($chrom2_name, $chrom2_start_bin, $chrom2_endin_bin) = ($chrom1_name, $chrom1_start_bin, $chrom1_endin_bin);
    }
    elsif ( scalar @chr_temp == 2 )
    {   
        ($chrom2_name, $chrom2_start_bin, $chrom2_endin_bin) = extract_location( $chr_temp[1] );
    }
    else
    {   
        print STDERR "Error parsing the chromosome bin information.\n\n";
        die get_use();
    }

    if ( !exists $chr_to_size{ $chrom1_name } )
    {
        die "Chromosome/Scaffold $chrom1_name not in the file.\n";
    }
    if ( !exists $chr_to_size{ $chrom2_name } )
    {
        die "Chromosome/Scaffold $chrom2_name not in the file.\n";
    }

    $chrom1_rowcol_num = int( $chr_to_size{ $chrom1_name } / $res ) + 1;
    $chrom2_rowcol_num = int( $chr_to_size{ $chrom2_name } / $res ) + 1;

    if ( $chrom1_start_bin == $UNINIT_VALUE ) { $chrom1_start_bin = 0; }
    if ( $chrom2_start_bin == $UNINIT_VALUE ) { $chrom2_start_bin = 0; }
    if ( $chrom1_endin_bin == $UNINIT_VALUE ) { $chrom1_endin_bin = $chrom1_rowcol_num - 1; }
    if ( $chrom2_endin_bin == $UNINIT_VALUE ) { $chrom2_endin_bin = $chrom2_rowcol_num - 1; }
}

if ($chrom1_name)
{
    print STDERR " Location1: $chrom1_name : $chrom1_start_bin, $chrom1_endin_bin\n";
}
if ($chrom2_name)
{
    print STDERR " Location2: $chrom2_name : $chrom2_start_bin, $chrom2_endin_bin\n";
}

my $inter_chrom_jump = 0;
if ( $inter_locn )
{
    my $chrom1_index = first { $sorted_chr_list[$_] eq $chrom1_name } 0 .. $#sorted_chr_list;
    my $chrom2_index = first { $sorted_chr_list[$_] eq $chrom2_name } 0 .. $#sorted_chr_list;
}
else
{
    if ( $chrom1_name && $chrom1_name ne $chrom2_name )
    {
        print STDERR "The file does not encode interchromosomal matrices.\n\n";
        exit 1;
    }
}

#Sanity checks
if ( $chrom1_start_bin < 0 || $chrom1_endin_bin < 0 || $chrom1_endin_bin < $chrom1_start_bin )
{
    die "The chromosome location/bin information is invalid.\n\n";
}

my @chrom1_list;
my @chrom2_list;

if ( !defined($location_str) && !defined($bin_str) )
{
    @chrom1_list = @sorted_chr_list;
    @chrom2_list = @sorted_chr_list;
}
else
{
    push @chrom1_list, $chrom1_name;
    push @chrom2_list, $chrom2_name;
}

foreach my $chrom1 ( @chrom1_list )
{
    foreach my $chrom2 ( @chrom2_list )
    {
        for (my $i = $chrom1_start_bin; $i <= $chrom1_endin_bin; $i++)
        {
            my @row_list;
            #intrachromosomal
            for (my $j = $chrom2_start_bin; $j <= $i; $j++)
            {
                push @row_list, get_values( $butlr_fp, $chr_to_locn{ $chrom1 }, $j, $i, $i );
            }
            if ($chrom2_endin_bin >= $i)
            {
                push @row_list, get_values( $butlr_fp, $chr_to_locn{ $chrom1 }, $i, $i, $chrom2_endin_bin );
            }
            #Print
            for (my $j = $chrom2_start_bin; $j <= $chrom2_endin_bin; $j++)
            {
                if (defined( $row_list[$j] ))
                {
                    print $row_list[$j];
                }
                if ($j < $chrom2_endin_bin)
                {
                    print "\t";
                }
                else
                {
                    print "\n";
                }
            }

        }
    }
}

close $butlr_fp;
