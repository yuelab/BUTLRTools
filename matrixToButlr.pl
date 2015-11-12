#! /usr/local/bin/perl -w

#Yanli Wang
#yanly.wang@gmail.com
#Version 1.2
#Code to convert Hi-C matrix (contact matrices in multiple files) into the Binary Upper TrianguLer MatRix (BUTLR) format

use strict;
use warnings;

#Find the path of the required modules
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/.");

#Required modules
use Getopt::Long qw(GetOptions);
use File::Basename;
use List::Util qw(any);
use Butlr;

my $version = "1.2";
my $assembly;
my $genome_size_filename;
my $matrix_list_filename;
my $resolution;
my $output_filename;
my $header_end_n = 1;
#my $i = 1; my $j = 2; my $v = 3; #Column # of i, j and Hi-C count value, respectively; may take up parameters in later verions

sub get_use()
{
    return  "BUTLR conversion tool version $version\n" .
            "Usage: perl $0 <REQUIRED> <OPTIONAL>\n\t" . 
            "-a <assembly of data> [REQUIRED]\n\t" .
            "\tPlease use assembly names recognized by UCSC Genome Browser.\n\t" . 
            "\tIf using new assembly, it is recommended to construct an assembly hub on UCSC and specify that URL, rather than default genome tracks.\n\t" .
            "-g <genome size file> [REQUIRED]\n\t" .
            "\tA file with chromosome/scaffold names and their sizes, delimited by tab, one per line, like the *.chrom.sizes files.\n\t" .
            "-m <file containing list of chr *tab* matrix file location> [REQUIRED]\n\t" . 
            "\tA file with chromosome/scaffold names and their corresponding matrix location/name, delimited by tab, one per line.\n\t" .
            "\tIf interchromosomal interactions are included, then the file could, for each line, list chrom1, chrom2 and file, deliminated by tab.\n\t" .
            "\tPlease make sure that the chromosome/scaffold here matches those from the genome size file.\n\t" . 
            "-r <resolution: bp as default, could enter Xk or Xm to specify resolutions in X kb or mb> [REQUIRED]\n\t" . 
            "-h <row number when header ends/matrix begins, 1-based; default: 1 (no header)> [OPTIONAL]\n\t" . 
            "-o <output butlr file> [OPTIONAL]\n";
}

GetOptions(
    'assembly|a=s' => \$assembly,
    'genome|g=s' => \$genome_size_filename,
    'matrix|m=s' => \$matrix_list_filename,
    'resolution|r=s' => \$resolution,
    'header_end|h=s' => \$header_end_n,
    'output|o=s' => \$output_filename,
) or die get_use();

unless (defined($assembly) && 
        defined($genome_size_filename) &&
        defined($matrix_list_filename) &&
        defined($resolution))
{
    die get_use();
}

my $res = get_resolution_in_bp($resolution);

#formulate the default output name if one is not provided by user
my ($basename, $dir, $ext) = fileparse($matrix_list_filename, qr/\.[^.]*/);
if (not defined $output_filename)
{
    $output_filename = "$dir$basename.$resolution" . ".butlr";
}

#read genome size file
if ( ! -e $genome_size_filename ) { die "$! ($genome_size_filename)\n" };
my %genome_size;
open(FILE, $genome_size_filename) or die "Error opening ($genome_size_filename)\n";
    while (my $line = <FILE>) 
    {
        chomp $line;
        my ($chr, $size) = split(/\t/, $line);
        $genome_size{$chr} = $size;
    }
close(FILE);

#sort by chromosome size, then by chromosome name (to put more "important" chromosomes first)
my @sorted_chr_list = ( sort { ($genome_size{$b} <=> $genome_size{$a}) || ($a cmp $b) } keys %genome_size );

my $interchrom_flag = 0;

#read matrix filenames
if ( ! -e $matrix_list_filename) { die "$! ($matrix_list_filename)\n" }
my %matrix_name;
open(FILE, $matrix_list_filename) or die "Error opening ($matrix_list_filename)\n";
    while (my $line = <FILE>)
    {
        chomp $line;
        my @recs = split(/\t/, $line);
        if ( scalar(@recs) == 2)
        {
            $matrix_name{ $recs[0] } = $recs[1];
        }
        elsif ( scalar(@recs) > 2 )
        {
            if ( $recs[0] eq $recs[1] )
            {
                $matrix_name{ $recs[0] } = $recs[2];
            }
            else
            {
                $interchrom_flag = 1;
                $matrix_name{ $recs[0] . "\t" . $recs[1] } = $recs[2];
            }
        }
        else
        {
            die "The matrix file does not seem to be in the correct format";
        }
    }
close(FILE);

open my $OUTFILE, ">", $output_filename or die "$! ($output_filename)\n";
binmode $OUTFILE;

#size of header (4-byte), accumulate size in this variable, then seek back to this location to overwite it once all headers are written
my $head_size = 0;
print $OUTFILE pack("L",  0);
$head_size += 4;

#version of the code that produced the format (fixed 16 1-byte)
my @clist = convert_from_char_to_bin($version);
while ( scalar(@clist) < 16 )
{
    push @clist, 0;  
}
print $OUTFILE pack("C*", @clist);
$head_size += 16;

#location of chromosome size/intrachromosomal
print $OUTFILE pack("L",  0);
$head_size += 4;

#location of interchromosomal information, 0 if does not exists
print $OUTFILE pack("L",  0);
$head_size += 4;

print STDERR "Assembly: " . $assembly . "\n";
#assembly
$assembly = trim($assembly);
print $OUTFILE pack("C*", convert_from_char_to_bin($assembly));
$head_size += length($assembly) + 1;

print STDERR "Resolution: " . $res . " bp\n";
#resolution (in bp)
print $OUTFILE pack("L",  $res);
$head_size += 4;

my $mcv = 0.0;
#Most common value in matrix (zero for now)
print $OUTFILE pack("f<",  $mcv);
$head_size += 4;

#Empty field holders x 4
print $OUTFILE pack("L",  0);
$head_size += 4;
print $OUTFILE pack("L",  0);
$head_size += 4;
print $OUTFILE pack("L",  0);
$head_size += 4;
print $OUTFILE pack("L",  0);
$head_size += 4;

#Denote the start of chromosome name, size, intrachromosomal information location 
seek ($OUTFILE, 20, 0); #0 = SEEK_SET
print $OUTFILE pack("L", $head_size);
seek ($OUTFILE,  0, 2); #2 = SEEK_END, go to end of file

#stores location of the first row for each intrachromosomal, first with 0, then come back to rewrite to these locations
my %intrachrom_to_row_header_location;
my %intrachrom_to_row_body_location;
foreach my $chr (@sorted_chr_list)
{
    if (exists $matrix_name{$chr})
    {
        if ( ! -e $matrix_name{$chr}) {unlink $output_filename; die "$! ($matrix_name{$chr})\n"}
        print $OUTFILE pack("C*", convert_from_char_to_bin($chr));
        $head_size += length($chr) + 1;
        print $OUTFILE pack("L",  $genome_size{$chr});
        $head_size += 4;
        $intrachrom_to_row_header_location{$chr} = $head_size;
        print $OUTFILE pack("Q",  0);
        $head_size += 8;
    }
}

seek ($OUTFILE, 0, 0);
print $OUTFILE pack("L", $head_size);
seek ($OUTFILE, 0, 2);

my $file_location = $head_size;
#
foreach my $chr (@sorted_chr_list)
{
    if (exists $matrix_name{$chr})
    {
        my $total_bin_num = int($genome_size{$chr} / $res) + 1;
        my $total_matrix  = ( $total_bin_num ** 2 + $total_bin_num ) / 2;
        print STDERR "Processing $chr : $total_bin_num length and $total_matrix bins with $matrix_name{$chr}\n";

        my @row_locations;
        if (open my $MXFILE, "<", $matrix_name{$chr})
        {
            my $n = 0;
            #Store sparse matrix in memory
            my $i = 0; #row number

            while (my $line = <$MXFILE>)
            {
                $n += 1;
                if ($n < $header_end_n)
                {
                    next;
                }
                if ($i >= $total_bin_num)
                {
                    warn "Larger matrix ($i >= $total_bin_num) than expected. Please make sure that the assembly genome file, matrix file and resolution are correct. Disregarding the excess rows.";
                    last;
                }
                chomp $line;
                $line = trim($line);
                my @recs = split(/\t/, $line);
                my $rec_size = scalar(@recs);

                if ( ($rec_size >= $total_bin_num) )
                {
                    @recs = splice(@recs, $rec_size - $total_bin_num, $total_bin_num);
                    #record row location
                    push @row_locations, $file_location;
                    #record col location
                    if ( any { $_ != $mcv } @recs )
                    {
                        for (my $j = $i; $j < $total_bin_num; $j++)
                        {
                            if ( $recs[$j] != $mcv )
                            {
                                print $OUTFILE pack("L", $j);
                                $file_location += 4;
                                print $OUTFILE pack("f<", $recs[$j]);
                                $file_location += 4;
                            }
                        }
                    }
                    #else
                #    {
                        #push @row_locations, 0;
                #    }
                }
                else { unlink $output_filename; die "Error: Incongruent matrix sizes\n"; }

                $i += 1;
            }
            close $MXFILE;
            push @row_locations, $file_location;
        }
        else {unlink $output_filename; die "$! ($matrix_name{$chr})\n";}
        if (scalar(@row_locations) < $total_bin_num)
        {
            unlink $output_filename; die "Error: Incongruent matrix sizes\n";
        }
        $intrachrom_to_row_body_location{ $chr } = $file_location;
        for (my $row = 0; $row < scalar( @row_locations ); $row++)
        {
            print $OUTFILE pack("Q", $row_locations[$row]);
            $file_location += 8;
        }
    } #if (exists $matrix_name{$chr})
} #foreach my $chr (@sorted_chr_list)

foreach my $chr (@sorted_chr_list)
{
    if (exists $matrix_name{$chr})
    {
        seek ($OUTFILE, $intrachrom_to_row_header_location{ $chr }, 0);
        print $OUTFILE pack("Q", $intrachrom_to_row_body_location{ $chr });
    }
}

close $OUTFILE;

print STDERR "Written $resolution" . " matrix as binary to $output_filename\n";
