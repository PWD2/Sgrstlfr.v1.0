#!/usr/bin/perl -w

use strict;
<<<<<<< HEAD
#use Data::Printer;

my $usage=<<USAGE;

    v1.0 2023.04 pangwending\@genomics.cn  draw dot graph from edge and nod file

=======
use Data::Printer;

my $usage=<<USAGE;

>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
    perl $0 <edge> <nod> <name> <outdir>

USAGE

my ($edg,$nod,$name,$outdir) = @ARGV;

die $usage if (@ARGV != 4);

open OUT,">$outdir/$name.1.dot" or die $!;
print OUT "graph $name\_dot_graph{\n";
print OUT '  rankdir=LR'."\n";
# print OUT '  dir=none'."\n";
my %edge_list;
my $genome_length=0;
open INA,$edg or die $!;
while(<INA>){
    chomp;
    next if (/^\#/);
    my ($line,$num) = split/\s+/,$_;
    next if ($num == 0);
    $genome_length +=$num;
    if ($num == 1){
        print OUT "  \"$line\:$num\"\[color=red\]\n";
    }
    $edge_list{$line} = $num;
}

my %line_hash;
my %cn_line;
open INB,$nod or die $!;
while(<INB>){
    chomp;
    my $anod = (split/\s+/,$_)[0];
    my ($ldex,$rdex) = split/\//,$anod;
    my $left_sight = judge_odd_even($ldex,'l');
    my $right_sight = judge_odd_even($rdex,'r');
    my $sort_left = re_sort_one_edge($left_sight);
    my $sort_right = re_sort_one_edge($right_sight);
    my $l_cn;my $r_cn;
    if (exists $edge_list{$sort_left}){
        $l_cn = $edge_list{$sort_left};
    }else{
        print "ERROR: The nod $anod left sight edge $sort_left is empty !!!!\n"
    }
    if (exists $edge_list{$sort_right}){
        $r_cn = $edge_list{$sort_right};
    }else{
        print "ERROR: The nod $anod right sight edge $sort_right is empty !!!!\n";
    }
    if ($l_cn and $r_cn){
        my $l_cn_line = "  \"$sort_left\:$l_cn\"--\"$sort_right\:$r_cn\"\[label=\"$anod\"\]";
        $line_hash{$l_cn_line} = '1';
    }
}

print OUT join("\n",keys %line_hash)."\n";
print OUT "  label = \"Chromosome length: $genome_length\"\n";
print OUT '}'."\n";
close OUT;
###################################
sub judge_odd_even{
    my ($input,$direct) = @_;
    my $yu = $input % 2;
    my $new_link;
    if ($yu == 1 ){
        my $asight = $input - 1;
        if ($direct eq 'l'){
            $new_link = $asight.'_'.$input;
        }else{
            $new_link = $input.'_'.$asight;
        }
    }else{
        my $asight = $input + 1;
        if ($direct eq 'r'){
            $new_link = $input.'_'.$asight;
        }else{
            $new_link = $asight.'_'.$input;
        }
    }
    return $new_link;
}
################################
sub re_sort_one_edge{
    my $input = $_[0];
    my $new_edge;
    my ($l,$r)=(split/_/,$input)[0,1];
    if ($l > $r){
        $new_edge = $r.'_'.$l;
        return $new_edge;
    }else{
        return $input;
    }
}

`dot -Tpdf $outdir/$name.1.dot -o $outdir/$name.1.dot_graph.pdf`;

