#!/usr/bin/perl -w

use strict;
use Data::Printer;
use Basicread qw/mersh_edge_region/;
use List::Util qw/pairs/;
use List::MoreUtils qw/uniq/;
use POSIX qw/floor/;
use Getopt::Long;

my ($edge,$loxpregion,$map) = @ARGV;
my ($samtools,$chrtype,$name,$outdir);
my $usage=<<USAGE;

Description: v1.0 2022.04.09 pwd 
             This script is purpose to get the breakpoint by split-reads information!!!
             (samtools path is must!)

Usage: perl $0 <edge> <loxpregion> <split_map> [-option]

            -t              :the chrmosome type [liner or cycle]
            -n              :the file name
            -s              :the path of samtools
            -o              :the output diectory

USAGE

die $usage if (!$edge||!$loxpregion||!$map);

GetOptions(
    "t:s" => \$chrtype,
    "n:s" => \$name,
    "s:s" => \$samtools,
    "o:s" => \$outdir
);

my %loxp_region = mersh_edge_region($edge,$loxpregion);
$samtools ||= '/hwfssz1/ST_BIOCHEM/P18Z10200N0255/PMO/F14ZQSYJSY1726/pangwending/some_tools/samtools-1.13/bin/samtools';
my $loxpsym = 'ATAACTTCGTATAATGTACATTATACGAAGTTAT';

###########Step1 make sure the range of the number ##################
my %new_loxp_region;
foreach my $keys (keys %loxp_region){
    my @set = @{$loxp_region{$keys}};
    next if ($set[0] == 0);
    my ($left,$right) = split/_/,$keys;
    my $mdian = $set[1]/2;
    my $plus;
    if ($mdian < 110){
        $plus = $mdian;
    }else{
        $plus = 110;
    }
    @{$new_loxp_region{$left}} = ($set[2],$set[2]+$plus - 1);
    @{$new_loxp_region{$right}} = ($set[3]-$plus + 1,$set[3]);   
}

############Step2 get the true start_point by loxpsym ###########################
open INA,"$samtools view $map |" or die $!;
my %pair_reads;
while(<INA>){
    chomp;
    my ($id,$flag,$start_point1,$map_quality,$my_seq) = (split/\s+/,$_)[0,1,3,4,9];
    # my @tranfer_flag = tentotwo($flag);
    # next if ($tranfer_flag[2] == 1);
    # next if ($map_quality < 30);
    my $check_start;
    if ($my_seq =~ /^$loxpsym/){
        $check_start = $start_point1 + 17;
    }elsif($my_seq =~ /$loxpsym$/){
        $check_start = $start_point1;
    }else{
        print "$id\n";
        print 'why not match the loxpsym???'."\n";
    }
    push @{$pair_reads{$id}},($check_start);
}
close INA;

###############Step3 locate the reads place in the edge range #################
my %new_pair_set;
foreach my $keys1 (keys %pair_reads){
    ##**repair
    my @reads_set = @{$pair_reads{$keys1}};
    #p(@reads_set);
    my $pair_num = scalar(@reads_set);
    if ($pair_num % 2 == 1 ){
        print 'extists possible bug!!!'."\n";
        p(@reads_set);
    }
    next if ($pair_num % 2 != 0);
    my $barcode_id;
    if ($keys1 =~ /\#(.+)/){
        $barcode_id = $1;
    }
    if ($pair_num == 2){
        my $a_l = corret_pe($reads_set[0],\%new_loxp_region);
        my $a_r = corret_pe($reads_set[1],\%new_loxp_region);
        my $this_nod = make_sure_pe($a_l,$a_r,$chrtype);
        push @{$new_pair_set{$barcode_id}},($this_nod) if ($this_nod ne 'undef');
    }else{
        foreach my $pair (pairs @reads_set){
            my ($key,$value) = @$pair;
            my $a_l = corret_pe($key,\%new_loxp_region);
            my $a_r = corret_pe($value,\%new_loxp_region);
            my $this_nod = make_sure_pe($a_l,$a_r,$chrtype);
            push @{$new_pair_set{$barcode_id}},($this_nod) if ($this_nod ne 'undef');
        }
    }
}

open OUT,">$outdir/$name.split.nod" or die $!;
foreach my $keys2 (keys %new_pair_set){
    my @no = @{$new_pair_set{$keys2}};
    print OUT $keys2."\t".join(",",@no)."\n";
}
close OUT;

##########################################
sub make_sure_pe{
    my ($a,$b,$type) = @_;
    my $a_nod = 'undef';
    if ($a ne 'undef' and $b ne 'undef'){
        if ($type eq 'liner'){
            if ($a != 0 and $b !=0){
                if ($a != $b){
                    $a_nod = $a.'/'.$b;
                }
            }
        }else{
            if ($a != $b){
                $a_nod = $a.'/'.$b;
            }
        }
    }
    return $a_nod;
}

############sub Ten to Two #################
sub tentotwo{
    my $num = $_[0];
    my @two;
    while($num >0){
        my $t =$num %2;
        push @two,$t;
        $num = floor($num/2);
        if ($num == 1){
            push @two,1;
            last;
        }else{
            next;
        }
    }
    return @two;
}

#######################################################################
sub corret_pe{
    my ($point,$range) = @_;
    my $num;
    foreach my $keys (keys %$range){
        my @range_set = @{${$range}{$keys}};
        if ( $range_set[0] <= $point and $point <= $range_set[1]){
            $num = $keys;
        }else{
            next;
        }
    }
    unless ($num){
        $num = 'undef';
    }
    return $num;
}





