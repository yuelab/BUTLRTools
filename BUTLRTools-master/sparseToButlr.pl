#! /usr/local/bin/perl -w

#Yanli Wang
#3dgenome.browser@gmail.com
#Version 1.3
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
use Butlr;

my $version = "1.3";
my $assembly;
my $genome_size_filename;
my $matrix_list_filename;
my $resolution;
my $output_filename;
my $header_end_n = 1;
my $i = 1; my $j = 2; my $v = 3; #Column # of i, j and Hi-C count value

sub get_use()
{
    return  "BUTLR conversion tool version $version\n" .
            "Converts contact matrices in coordinate list to BUTLR\n\n" .
            "Usage: perl $0 <REQUIRED> <OPTIONAL>\n\t" . 
            "-a <assembly of data> [REQUIRED]\n\t" .
            "\tPlease use assembly names recognized by UCSC Genome Browser.\n\t" . 
            "\tIf using new assembly, it is recommended to construct an assembly hub on UCSC and specify that URL, rather than default genome tracks.\n\n\t" .
            "-g <genome size file> [REQUIRED]\n\t" .
            "\tA file with chromosome/scaffold names and their sizes, delimited by tab, one per line, like the *.chrom.sizes files.\n\n\t" .
            "-m <file containing list of chr *tab* matrix file location> [REQUIRED]\n\t" . 
            "\tA file with chromosome/scaffold names and their corresponding matrix location/name, delimited by tab, one per line.\n\t" .
            "\tIf interchromosomal interactions are included, then the file could, for each line, list chrom1, chrom2 and file, deliminated by tab.\n\t" .
            "\tPlease make sure that the chromosome/scaffold here matches those from the genome size file.\n\n\t" . 
            "-r <resolution: bp as default, could enter Xk or Xm to specify resolutions in X kb or mb> [REQUIRED]\n\n\t" . 
            "-i <column # of row/chrom1 information, 1-based; default: 1> [OPTIONAL]\n\n\t" .
            "-j <column # of column/chrom2 information, 1-based; default: 2> [OPTIONAL]\n\n\t" .
            "-v <column # of Hi-C value information, 1-based; default: 3> [OPTIONAL]\n\n\t" .
            "-h <row number when header ends/matrix begins, 1-based; default: 1 (no header)> [OPTIONAL]\n\n\t" . 
            "-o <output butlr file> [OPTIONAL]\n\n\n";
}

GetOptions(
    'assembly|a=s' => \$assembly,
    'genome|g=s' => \$genome_size_filename,
    'matrix|m=s' => \$matrix_list_filename,
    'resolution|r=s' => \$resolution,
    'i=s' => \$i,
    'j=s' => \$j,
    'v=s' => \$v,
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
    $output_filename = "$dir$basename.$resolution" . ".btr";
}

#read genome size file
if ( ! -e $genome_size_filename ) { die "$! ($genome_size_filename)\n" };
my %genome_size;
open(FILE, $genome_size_filename) or die "Error opening ($genome_size_filename)\n";
    while (my $line = <FILE>) 
    {
        chomp $line;
        my ($chr, $size) = split(' ', $line);
        $genome_size{$chr} = $size;
    }
close(FILE);

#sort by chromosome size, then by chromosome name (to put more "important" chromosomes first)
############################################################################################################
my @sorted_chr_list = ( sort { ($genome_size{$b} <=> $genome_size{$a}) || ($a cmp $b) } keys %genome_size );
############################################################################################################

my $interchrom_flag = 0;

#read matrix filenames
if ( ! -e $matrix_list_filename) { die "$! ($matrix_list_filename)\n" }
my %matrix_name;
open(FILE, $matrix_list_filename) or die "Error opening ($matrix_list_filename)\n";
    while (my $line = <FILE>)
    {
        chomp $line;
        my @recs = split(' ', $line);
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

my %interchrom_to_row_header_location;
my %interchrom_to_row_body_location;
if ( $interchrom_flag )
{
my $interchrom_loctn = $head_size;

foreach my $c1 (0 .. $#sorted_chr_list)
{
    foreach my $c2 (0 .. $c1)
    {
        my $chrom1 = $sorted_chr_list[$c1];
        my $chrom2 = $sorted_chr_list[$c2];
        if ( $chrom1 eq $chrom2 ) { next; }

        if ( exists $matrix_name{ $chrom1 . "\t" . $chrom2 } || exists $matrix_name{ $chrom2 . "\t" . $chrom1 } )
        {
            my $fname;
            if ( exists $matrix_name{ $chrom1 . "\t" . $chrom2 } )
            {
                $fname = $matrix_name{ $chrom1 . "\t" . $chrom2 };
            }
            if ( exists $matrix_name{ $chrom2 . "\t" . $chrom1 } )
            {
                $fname = $matrix_name{ $chrom2 . "\t" . $chrom1 };
            }

            if ($fname)
            {
                my $key;
                if ( is_chrom1_ahead(\%genome_size, $chrom1, $chrom2) )
                {
                    $key = $chrom2 . "\t" . $chrom1;
                }
                else
                {
                    $key = $chrom1 . "\t" . $chrom2;
                }

                if ( ! -e $fname ) {unlink $output_filename; die "$! ($fname)\n"}
                print $OUTFILE pack("C*", convert_from_char_to_bin($key));
                $head_size += length($key) + 1;
                $interchrom_to_row_header_location{ $key } = $head_size;
                print $OUTFILE pack("Q",  0); #fill this later
                $head_size += 8;
            }
        }
    }
}

seek($OUTFILE, 24, 0);
print $OUTFILE pack("L", $interchrom_loctn);
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

        my @sparse_list = ();
        if (open my $MXFILE, "<", $matrix_name{$chr})
        {
            my $n = 0;
            while (my $line = <$MXFILE>)
            {
                $n += 1;
                if ($n < $header_end_n)
                {
                    next;
                }
                chomp $line;
                $line = trim($line);
                my @recs = split(' ', $line);

                my $a = $recs[$i-1] / $res;
                my $b = $recs[$j-1] / $res;
                if ( $a > $b ) #Makes sure that values under index i is always <= j
                {
                    $a = $recs[$j-1] / $res;
                    $b = $recs[$i-1] / $res;
                }
                push( @{$sparse_list[$n-$header_end_n]}, ($a, $b, $recs[$v-1]) );
            }
            close $MXFILE;
        }
        else {unlink $output_filename; die "$! ($matrix_name{$chr})\n";} 

        #Sort by i, then j
        @sparse_list = sort { ( $a->[0] <=> $b->[0] ) || ( $a->[1] <=> $b->[1] ) } @sparse_list;
        my @row_locations;
        my $x = 0;
        for (my $row = 0; $row < $total_bin_num; $row++)
        {
            if (defined $sparse_list[$x])
            {
                if ($sparse_list[$x][0] != $row)
                {
                    push @row_locations, 0;
                }
                else
                {
                    push @row_locations, $file_location;
                }
                while (defined $sparse_list[$x] and $sparse_list[$x][0] == $row)
                {
                    print $OUTFILE pack("L",  $sparse_list[$x][1]);
                    $file_location += 4;
                    print $OUTFILE pack("f<", $sparse_list[$x][2]);
                    $file_location += 4;
                    $x += 1;
                }
            }
            else
            {
                push @row_locations, 0;
            }
        }
        push @row_locations, $file_location;
        my $row_number = scalar(@row_locations);
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

seek ($OUTFILE, 0, 2);

if ( $interchrom_flag )
{
foreach my $c1 ( 0 .. $#sorted_chr_list )
{
    foreach my $c2 ( 0 .. $c1 )
    {
        my $chrom1 = $sorted_chr_list[$c1];
        my $chrom2 = $sorted_chr_list[$c2];

        #Skip intrachromosomal datasets
        if ( $chrom1 eq $chrom2 ) { next; }

        my $rev_flag = 0;
        my $key;
        if ( exists $matrix_name{  $chrom2 . "\t" . $chrom1 } )
        {
            $key = $chrom2 . "\t" . $chrom1;
            $rev_flag = 1;
        }
        else
        {
            $key = $chrom1 . "\t" . $chrom2;
            $rev_flag = 0;
        }

        if ( exists $matrix_name{ $key } )
        {
            my $fname = $matrix_name{ $key };

            if ($fname)
            {
                print STDERR $key . ": ";

                my $c_chrom = $chrom2;
                my $r_chrom = $chrom1;               
                if ( is_chrom1_ahead(\%genome_size, $chrom1, $chrom2) )
                {
                    my $c_chrom = $chrom1;
                    my $r_chrom = $chrom2;
                }

                my $chrom_a_bin = int($genome_size{ $c_chrom } / $res) + 1;
                my $chrom_b_bin = int($genome_size{ $r_chrom } / $res) + 1;
                my $total_bin   = $chrom_a_bin * $chrom_b_bin;

                print STDERR "Processing interchromosomal $r_chrom ($chrom_b_bin rows) x $c_chrom ($chrom_a_bin columns) [total: $total_bin bins] with $fname\n";


                my @sparse_list = ();
                if (open my $MXFILE, "<", $fname)
                {
                    my $n = 0;
                    while (my $line = <$MXFILE>)
                    {
                        $n += 1;
                        if ($n < $header_end_n)
                        {
                            next;
                        }
                        chomp $line;
                        $line = trim($line);
                        my @recs = split(' ', $line);

                        if ( $rev_flag )
                        {
                            if ( $recs[$i-1] > $genome_size{$c_chrom} )
                            {
                                unlink $output_filename;
                                die "The designated bin is larger than chromosome size. Please make sure that the order of chromosomes listed in the matrix list file is the same as the coordinated list file 1\n";
                            }
                            if ( $recs[$j-1] > $genome_size{$r_chrom} )
                            {
                                unlink $output_filename;
                                die "The designated bin is larger than chromosome size. Please make sure that the order of chromosomes listed in the matrix list file is the same as the coordinated list file 2\n";
                            }
                        }
                        else
                        {
                            if ( $recs[$j-1] > $genome_size{$c_chrom} )
                            {
                                unlink $output_filename;
                                die "The designated bin is larger than chromosome size. Please make sure that the order of chromosomes listed in the matrix list file is the same as the coordinated list file 1\n";
                            }
                            if ( $recs[$i-1] > $genome_size{$r_chrom} )
                            {
                                unlink $output_filename;
                                die "The designated bin is larger than chromosome size. Please make sure that the order of chromosomes listed in the matrix list file is the same as the coordinated list file 2\n";
                            }
                        }

                        my $a = $recs[$i-1] / $res;
                        my $b = $recs[$j-1] / $res;
                        if ( $rev_flag )
                        {
                            $a = $recs[$j-1] / $res;
                            $b = $recs[$i-1] / $res;
                        }
                        push( @{$sparse_list[$n-$header_end_n]}, ($a, $b, $recs[$v-1]) );
                    }
                     close $MXFILE;
                }
                else {unlink $output_filename; die "$! ($fname)\n";}

		#Sort by i, then j
                @sparse_list = sort { ( $a->[0] <=> $b->[0] ) || ( $a->[1] <=> $b->[1] ) } @sparse_list;
                my @row_locations;
                my $x = 0;
                for (my $row = 0; $row < $chrom_a_bin; $row++)
                {
                    if (defined $sparse_list[$x])
                    {
                        if ($sparse_list[$x][0] != $row)
                        {
                            push @row_locations, 0;
                        }
                        else
                        {
                            push @row_locations, $file_location;
                        }
                        while (defined $sparse_list[$x] and $sparse_list[$x][0] == $row)
                        {
                            print $OUTFILE pack("L",  $sparse_list[$x][1]);
                            $file_location += 4;
                            print $OUTFILE pack("f<", $sparse_list[$x][2]);
                            $file_location += 4;
                            $x += 1;
                        }
                    }
                    else
                    {
                        push @row_locations, 0;
                    }
                }
                push @row_locations, $file_location;
                if ( is_chrom1_ahead(\%genome_size, $chrom1, $chrom2) )
                {
                    $interchrom_to_row_body_location{ $chrom2 . "\t" . $chrom1 } = $file_location;
                }
                else
                {
                    $interchrom_to_row_body_location{ $chrom1 . "\t" . $chrom2 } = $file_location;
                }

                for (my $row = 0; $row < scalar( @row_locations ); $row++)
                {
                    print $OUTFILE pack("Q", $row_locations[$row]);
                    $file_location += 8;
                }
            }
        }
    } #foreach my $c2 ($c1 .. $#sorted_chr_list)
} #foreach my $c1 (0 .. $#sorted_chr_list)
} #if ( $interchrom_flag )


foreach my $c1 ( 0 .. $#sorted_chr_list )
{
    foreach my $c2 ( 0  .. $c1 )
    {
        my $chrom1 = $sorted_chr_list[$c1];
        my $chrom2 = $sorted_chr_list[$c2];
        if ( $chrom1 eq $chrom2 ) { next; }

        my $key;
        if ( is_chrom1_ahead(\%genome_size, $chrom1, $chrom2) )
        {
            $key = $chrom2 . "\t" . $chrom1;
        }
        else
        {
            $key = $chrom1 . "\t" . $chrom2;
        }

        if ( (exists $matrix_name{  $chrom1 . "\t" . $chrom2 } && -e $matrix_name{  $chrom1 . "\t" . $chrom2 }) ||
             (exists $matrix_name{  $chrom2 . "\t" . $chrom1 } && -e $matrix_name{  $chrom2 . "\t" . $chrom1 }) )
        {   
            my $fname = $matrix_name{ $key };
            seek ($OUTFILE, $interchrom_to_row_header_location{ $key }, 0);
            print $OUTFILE pack("Q", $interchrom_to_row_body_location{ $key });
        }
    }
}

close $OUTFILE;

print STDERR "Written $resolution" . " matrix as binary to $output_filename\n";
