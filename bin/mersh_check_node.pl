#!/usr/bin/perl -w

use strict;
use Data::Printer;
use List::MoreUtils qw/uniq/;
use Cwd qw/getcwd/;

#This script is order to mersh the two node file and re-check the node !

my ($o_node,$pe_node,$split_node,$name,$outdir) = @ARGV;

my $usage=<<USAGE;

Usage: perl $0 <o_node> <pe_node> <split_node> <name> <outdir>

        o_node              :the corret node set
        pe_node             :the pe node set
        split_node          :the split reads node set
        name                :output file name [default:test]
        outdir              :output file diectiony [outdir]

USAGE

die $usage if (!$o_node ||!$pe_node||!$split_node);

$name ||='test';
$outdir ||=getcwd;

my @my_corret_set;
open INA,$o_node or die $!;
while(<INA>){
    chomp;
    my $node = (split/\s+/,$_)[0];
    push @my_corret_set,$node;
}

close INA;

my %sum_set;
open INB,$pe_node or die $!;
while(<INB>){
    chomp;
    my ($barcode,$ano) = split/\s+/,$_;
    my @ano_set = split/\,/,$ano;
    push @{$sum_set{$barcode}},(@ano_set);
}
close INB;

open INC,$split_node or die $!;
while(<INC>){
    chomp;
    my ($barcode_id,$aano) = split/\s+/,$_;
    my @aano_set = split/\,/,$aano;
    push @{$sum_set{$barcode_id}},(@aano_set);
}

close INC;

open OUT,">$outdir/$name.mersh.nod" or die $!;
foreach my $keys (keys %sum_set){
    my @this_nod = @{$sum_set{$keys}};
    my @a_new_set;
    foreach my $keys1 (@this_nod){
        my ($left,$right) = split/\//,$keys1;
        my $new_node;
        if ($left > $right){
            $new_node = $right.'/'.$left;
        }else{
            $new_node = $keys1
        }
        if (grep {$new_node eq $_} @my_corret_set){
            push @a_new_set,$new_node;
        }else{
            print '#exsits erro node left!'."\n";
            print '#you should repair the pe-nod and split-nod or ignore !'."\n";
            print $new_node."\n";
        }
    }
    my @uniq_set = uniq(@a_new_set);
    if($uniq_set[0]){
        print OUT "$keys\t".join(",",@uniq_set)."\n";
    }
}

close OUT;
