#Yanli Wang
#3dgenome.browser@gmail.com
#Version 1.3

package Butlr;
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw( convert_from_char_to_bin convert_from_bin_to_char get_resolution_in_bp trim is_chrom1_ahead );

# these are exported by default.
our @EXPORT    = qw( convert_from_char_to_bin convert_from_bin_to_char get_resolution_in_bp trim is_chrom1_ahead );

#Parameter: String
#Return: ArrayList of ASCII values of each character in the string and 0 (designated null terminator)
sub convert_from_char_to_bin
{
    my @ascii_list = ();
    my @char_list = split("", $_[0]);
    for (my $i = 0; $i < scalar(@char_list); $i++)
    {
        push @ascii_list, ord $char_list[$i];
    }
    push @ascii_list, 0;
    return @ascii_list;
}

sub convert_from_bin_to_char
{
    my @ascii_list = @{$_[0]};
    my $string = '';
    for (my $i = 0; $i < scalar(@ascii_list); $i++)
    {
        if ($ascii_list[$i] == 0)
        {
            last;
        }
        $string .= chr $ascii_list[$i];
    }
    return $string;
}

sub get_resolution_in_bp
{
    my $s = $_[0];
    my $res = $s;
    my $i = index(lc($s), 'k');
    if ($i >= 0)
    {
        $res = substr($s, 0, $i) * 1000;
    }
    else
    {
        $i = index(lc($s), 'm');
        if ($i >= 0)
        {
            $res = substr($s, 0, $i) * 1000000;
        }   
    }
    return $res;
}

#http://perlmaven.com/trim
sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

#Returns true if chrom1 is prioritized over chrom2 (in order to reduce redundancy); else if otherwise
#Input:
#   1: hashtable with key: name of chromosome/scaffold and value: size
#   2: name of chrom1
#   3: name of chrom2
#Order of priority determined by chromosome size (larger chrom with more prioritiy) and then (if tie) alphabetically sorted name of chromosome/name
sub is_chrom1_ahead
{
    my %genome_size = %{$_[0]};
    my $chrom1 = $_[1];
    my $chrom2 = $_[2];

    if ( $chrom1 eq $chrom2 )
    {
        return 0;
    }

    if ( $genome_size{$chrom1} > $genome_size{$chrom2} )
    {
        return 1;
    }

    if ( $genome_size{$chrom1} == $genome_size{$chrom2} )
    {
        if ( ($chrom1 cmp $chrom2) < 0)
        {
            return 1;
        }
    }
    return 0;
}
1;
