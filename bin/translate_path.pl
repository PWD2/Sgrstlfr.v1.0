#!/usr/bin/perl -w

use strict;
use Basicread qw/mersh_edge_region read_fasta/;

my ($path,$region_file,$edge_file,$fasta) = @ARGV;
my ($chr_id,$name,$outdir);

my $usage=<<USAGE;

    Description: This script is translate the edge-nod path to fasta

    Usage: perl $0 <path> <region> <edge> <fasta> [-option]
            
            <path>      :the path which want to restructure [list or alon is ok]
            <region>    :this syn-chrmosome's region file
            <edge>      :this syn-chromosome's edge file
            <fasta>     :the fasta file
            -c          :the chrmosome name [must]
            -n          :the file name
            -o          :the output directory

USAGE

GetOptions(
    "c:s" =>\$chr_id,
    "n:s" =>\$name,
    "o:s" => \$outdir
);

die $usage if (@ARGV !=4 ||!$chr_id);
my @path_set;
if (-e $path){
    open INA,$path or die $!;
    while(<INA>){
        chomp;
        push @path_set,$_;
    }
}else{
    push @path_set,$path;
}
my $num = 1;
foreach my $solution (@path_set){
    open OUTA,">$outdir/$name.$num.fa" or die $!;
    my $this_fasta = translate_path($edge_file,$region_file,$fasta,$solution,$chr_id);
    my $length = length($this_fasta);
    print OUTA ">$solution#$length"."\n";
    for (my $i=0;$i<$length;$i=$i+60){
        my $line_fa = substr($this_fasta,$i,60);
        print OUTA $line_fa."\n";
    }
    close OUTA;
    $num++;
}

#################################################################
sub translate_path{
    my ($edg,$region_set,$fasta,$path,$chrid) = @_;
    my %edge_region_set = mersh_edge_region($edg,$region_set);
    my $seq = read_fasta($fasta,$chrid);
    my $new_seq='';
    my @path_edge = split/\//,$path;
    foreach my $a_edge (@path_edge){
        my @edge_index = split/_/,$a_edge;
        if ($edge_index[0] > $edge_index[1]){
            my $new_edge = $edge_index[1].'_'.$edge_index[0];
            my $start_point = ${$edge_region_set{$new_edge}}[2];
            my $length = ${$edge_region_set{$new_edge}}[1];
            my $con_seq = substr($seq,$start_point-1,$length);
            my $rev_seq = reverse($con_seq);
            $rev_seq =~ tr/ATCG/TAGC/;
            $new_seq = $new_seq.$rev_seq;
        }else{
            my $start_point = ${$edge_region_set{$a_edge}}[2];
            my $length = ${$edge_region_set{$a_edge}}[1];
            my $con_seq = substr($seq,$start_point-1,$length);
            $new_seq = $new_seq.$con_seq;
        }
    }
    return $new_seq;
}

