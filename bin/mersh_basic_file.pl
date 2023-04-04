#!/usr/bin/perl -w

use strict;
<<<<<<< HEAD
#use Data::Printer;
=======
use Data::Printer;
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
use Basicread qw(read_edge read_region read_depth);
use List::Util qw/max min sum/;
use Cwd qw/getcwd/;

my $usage=<<USAGE;

Desciption: This script is used to mersh the loxp_region、edge、and depthsingle !

Usage: perl $0 <edge> <region> <depthsite> <chrid> <name> <outdir>

USAGE

my ($edge,$region,$depth,$chr,$name,$outdir)=@ARGV;

die $usage if (!$edge||!$region||!$depth||!$chr);

print 'Start time: '.localtime()."\n";
$name ||='test';
$outdir ||=getcwd;

my %edge_hash = read_edge($edge);
my %region_hash = read_region($region);
my @depth_set = read_depth($depth,$chr);

my %result;
foreach my $keys (keys %edge_hash){
    my @aset = @{$edge_hash{$keys}};
    #next if ($aset[0] == 0);
    my @index_set =split/\,/,$aset[1]; 
    if (@index_set == 1){
        my @region_range = @{$region_hash{$index_set[0]}};
        my $length = $region_range[2] - $region_range[1] +1;
        my $cover_depth;
        if ($aset[0] == 0){
            $cover_depth = '0.00';
        }else{
            my @the_depth_range = @depth_set[$region_range[1]..$region_range[2]-1];
            my $sum_depth = sum (@the_depth_range);
            $cover_depth = sprintf("%.2f",$sum_depth/$length);
        }
        @{$result{$keys}} =($aset[0],$length,$cover_depth,$region_range[1],$region_range[2],$aset[1]); 

    }else{
        my $max_index = max(@index_set);
        my $min_index = min(@index_set);
        my $max_point = ${$region_hash{$max_index}}[2];
        my $min_point = ${$region_hash{$min_index}}[1];
        my $length = $max_point - $min_point +1 ;
        my $cover_depth;
        if ($aset[0] == 0){
            $cover_depth = '0.00';
        }else{
            my @the_depth_range = @depth_set[$min_point..$max_point-1];
            my $sum_depth = sum (@the_depth_range);
            $cover_depth = sprintf("%.2f",$sum_depth/$length);
        }
        @{$result{$keys}} =($aset[0],$length,$cover_depth,$min_point,$max_point,$aset[1]);
    }
}

open OUTA,">$outdir/$name.total" or die $!;
print OUTA "#edge\t#CN\t#length\t#cover_depth\t#start_point\t#end_point\t#loxp_index\n";
my $total_length = 0 ;
my $total_cn = 0;
my $uniq_edge = 0;
my $dele_edge = 0;
my $edge_type = keys %result;

foreach my $keys (sort {(split/_/,$a)[0] <=> (split/_/,$b)[0]} keys %result){
    my @set = @{$result{$keys}};
    my $true_length = $set[0]*$set[1];
    $total_length = $total_length + $true_length;
    $total_cn = $total_cn + $set[0];
    if ($set[0] == 0){
        $dele_edge++;
    }elsif ($set[0] == 1){
        $uniq_edge++;
    }
    print OUTA "$keys\t".join("\t",@set)."\n";
}
print OUTA "\n";
print OUTA '#The SCRaMbLE genome length:'."$total_length\n";
print OUTA '#The total CN number:'."$total_cn\n";
print OUTA '#The edge type number:'."$edge_type\n";
print OUTA '#The uniq edge number:'."$uniq_edge\n";
print OUTA '#The delete edge number:'."$dele_edge\n";

close OUTA;
print 'End time: '.localtime()."\n";
print '###This work is finished!'."\n";







