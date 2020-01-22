#!/usr/bin/perl -w

# reformat LAT run processing times

# input format
# MET(seconds) SLAC(hours) NASA(hours)

# usage: 
# ./dates.pl < allruns.YYYYMMDD > datetimes.YYYYMMDD

#use POSIX qw(strftime);
use Time::Local;

# offet between Fermi MET seconds and Unix seconds (neglecting leap seconds)
$m2001 = 978307200;
$junk = 0;

print "     YMD      HMS      FYear       DOY        Unix      MET    SLAC   NASA   TOTAL\n"; 
while (<>) {
    ($met,$slac,$nasa) = split;
    $met = int($met + 0.5);   # do this because Perl "int" truncates towards 0 instead of rounding
    $ttot = $slac + $nasa;
    $unix = $met + $m2001;
    ($s,$m,$h,$mday,$mon,$year,$junk,$doy) = gmtime($unix);
    $year += 1900;
    $mon++;
    $doy += $h/24.0 + $m/1440.0 + $s/86400.0;
    $ymd = sprintf "%4d-%02d-%02d",$year,$mon,$mday;
    $hms = sprintf "%02d:%02d:%02d",$h,$m,$s;
    $ylen = ($year % 4)? 365.0 : 366.0;
    $fy = $year + $doy/$ylen;
#    printf "%.3f %.3f = $y/$doy $h:$m:$s\n", $fy,$slac;
    printf "%s %s %.6f %8.4f %s %s %.4f %.4f %.4f\n", $ymd,$hms,$fy,$doy,$unix,$met,$slac,$nasa,$ttot;
}

