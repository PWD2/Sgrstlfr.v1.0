#!1/usr/bin/perl -w

use strict;
use Basicread qw/mersh_edge_region/;
use List::MoreUtils qw/uniq/;
use Data::Printer;
use Getopt::Long;
use Cwd qw/getcwd/;

my ($edge,$region,$synmap) = @ARGV;
my ($chr_type,$insert_size,$overlap,$name,$outdir);

my $usage=<<USAGE;

Description: v1.0 2022.04.07 pwd 
             This script is purpose to stat the breakpoint by pair-end !

Usage:perl $0 <edge> <region> <synmap>

        -n          : the name of the output file [defacult:test]
        -ct         : the chr type of input file [liner or cycle]
        -isz        : the insert_size of the pair-end seq [defacult:500]
        -overlap    : the overlap of reads in range will be left [defacult:0.8]
        -o          : the output directory

USAGE

die $usage if (!$edge||!$region||!$synmap);

GetOptions(
    "n:s" => \$name,
    "ct:s" => \$chr_type,
    "isz:f" => \$insert_size,
    "overlap:f" =>\$overlap,
    "o:s" => \$outdir
);

$name ||='test';
$overlap ||= 0.8;
$insert_size ||= 500;
$outdir ||=getcwd;

my %edge_region = mersh_edge_region($edge,$region);
my %range_hash;
my $index = 0;
foreach my $keys (keys %edge_region){
    my @set = @{$edge_region{$keys}};
    next if ($set[0] == 0);
    next if ($set[1] < $insert_size); ### the erro range maybe produce bug !!!
    my ($left,$right) = split/_/,$keys;
    my $median = $set[1]/2;
    my $plus;
    if ($insert_size > $median){
        $plus = $median;
    }else{
        $plus = $insert_size;
    }
    push @{$range_hash{$left}},($set[2],$set[2] + $plus -1) ;   #1 3 5 7 9 
    push @{$range_hash{$right}},($set[3] - $plus + 1,$set[3]) ;    #2 4 6 8 10
    $index++;
}
#p(%range_hash);exit;

open INA,"gzip -dc $synmap |" or die $!;
my %new_pair_set;
while(<INA>){
    chomp;
    my ($id,$start_point1,$end_point1,$start_point2,$end_point2) = (split/\s+/,$_)[0,2,3,6,7];
    my $barcode;
    next if ($start_point2 - $start_point1 == 0);
    if ($id =~ /\#(.+)/){
        $barcode = $1;
    }
    my $id_l = pe_belong_range($start_point1,$end_point1,$overlap,\%range_hash);
    my $id_r = pe_belong_range($start_point2,$end_point2,$overlap,\%range_hash);
    # print 'l:'.$id_l."\t".'r:'.$id_r."\n";
    if ($id_l ne 'undef' and $id_r ne 'undef'){
        if ($chr_type eq 'liner'){
            if ($id_l != 0 and $id_r !=0){
                if ($id_l != $id_r){
                    my $a_nod = $id_l.'/'.$id_r;
                    push @{$new_pair_set{$barcode}},$a_nod;
                }
            }
        }elsif ($chr_type eq 'cycle'){
                if ($id_l != $id_r){
                    my $a_nod = $id_l.'/'.$id_r;
                    push @{$new_pair_set{$barcode}},$a_nod;
                }
        }else{
            print "Please input the corret chr-type !!!\n";
            exit;
        }
    }
}

open OUT,">$outdir/$name.pe.nod" or die $!;
foreach my $keys1 (keys %new_pair_set){
    my @this_set = @{$new_pair_set{$keys1}};
    my @uniq_set = uniq (@this_set);
    print OUT "$keys1\t".join(",",@uniq_set)."\n";
}
close OUT;

####################################################
sub pe_belong_range{
    my ($read_a,$read_b,$include,$a_hash) = @_;
    my $result_line;
    foreach my $line (keys %$a_hash){
        my $range_a = ${${$a_hash}{$line}}[0];
        my $range_b = ${${$a_hash}{$line}}[1];
        if ($range_a <= $read_a and $read_b <= $range_b){
            $result_line = $line;
        }elsif($read_a > $range_a and $range_b > $read_a){
            my $precent = ($range_b - $read_a + 1)/100;
            if ($precent >= $include){
                $result_line = $line;
            }
        }elsif($read_b < $range_b and $range_a < $read_b){
            my $precent = ($read_b - $range_a + 1)/100;
            if ($precent >= $include){
                $result_line = $line;
            }
        }else{
            next;
        }
    }
    if ($result_line){
        return $result_line;
    }else{
        return 'undef';
    }
}

