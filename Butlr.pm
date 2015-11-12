package Butlr;
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw( convert_from_char_to_bin convert_from_bin_to_char get_resolution_in_bp trim );

# these are exported by default.
our @EXPORT    = qw( convert_from_char_to_bin convert_from_bin_to_char get_resolution_in_bp trim );

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

1;
