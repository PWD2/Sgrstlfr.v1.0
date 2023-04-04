#!/usr/bin/perl -w

use strict;
use List::Util qw/max/;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;

    Descript: This script is fix the cnv number by the average depth and nod file !!!
              It can be reduce the conflict between cnv and breakpoint, but you should 
              keep the recombination loxPsym breakpoint is accurate enough.

    Usage:  perl $0 <edge> <nod> [-option]

            -mdep           : the average seqence depth in all genome(or you set in cnv estimation)
            -n              : the output file name [default: test]
            -o              : the output directory [default: ./]

USAGE


my ($edge_list,$nod_list) = @ARGV;
my ($aver_depth,$name,$outdir);

GetOptions(
    "mdep:s" => \$aver_depth,
    "n:s" => \$name,
    "o:s" => \$outdir
);

$name ||='test';
$outdir ||=getcwd;

die $usage if (!$edge_list||!$nod_list||!$aver_depth);

my %nod_index_hash;
open INB,$nod_list or die $!;
while(<INB>){
    chomp;
    my $the_nod = (split/\s+/,$_)[0];
    my ($ldex,$rdex) = split/\//,$the_nod;
    push @{$nod_index_hash{$ldex}},($the_nod);
    push @{$nod_index_hash{$rdex}},($the_nod);
}
close INB;

open OUT,">$outdir/$name.edg2" or die $!;
open INA,$edge_list or die $!;
while(<INA>){
    chomp;
    my @this_set = split/\s+/,$_;
    if ($this_set[1] == 0 and $this_set[3] < $aver_depth){
        my ($left,$right) = split/_/,$this_set[0];
        my ($lnum,$rnum);
        if (exists $nod_index_hash{$left}){
            $lnum = scalar(@{$nod_index_hash{$left}});
        }else{
            $lnum = 0;
        }
        if (exists $nod_index_hash{$right}){
            $rnum = scalar(@{$nod_index_hash{$right}});
        }else{
            $rnum = 0;
        }
        if ($lnum == 0 and $rnum == 0){
            print OUT $_."\n";
        }else{
            my $max = max($lnum,$rnum);
            print OUT $this_set[0]."\t".$max."\t".join("\t",@this_set[2..4])."\n";
        }
    }else{
        print OUT $_."\n";
    }
}
close INA;
close OUT;
# `mv "$outdir/$name.edg2" "$outdir/$name.edg"`;


