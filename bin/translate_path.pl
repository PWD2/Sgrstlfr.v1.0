#!/usr/bin/perl -w

use strict;
use Basicread qw/mersh_edge_region read_fasta/;
use Getopt::Long;
use List::MoreUtils qw/indexes/;
use Cwd;

my ($path,$region_file,$edge_file,$fasta) = @ARGV;
my ($chr_id,$chr_type,$r_type,$name,$outdir);

my $usage=<<USAGE;

    Description: This script is translate the edge-nod path to fasta

    Usage: perl $0 <path> <region> <edge> <fasta> [-option]
            
            <path>      :the path which want to restructure [list or alon is ok]
            <region>    :this syn-chrmosome's region file
            <edge>      :this syn-chromosome's edge file
            <fasta>     :the fasta file
            -a          :[all | first]produece the all fasta file or first one [default:first]
            -c          :the chromosome name [must]
            -t          :[linear | cycle ] the chromosome type 
            -n          :the file name
            -o          :the output directory

USAGE

GetOptions(
    "a:s" => \$r_type,
    "c:s" =>\$chr_id,
    "t:s" => \$chr_type,
    "n:s" =>\$name,
    "o:s" => \$outdir
);

die $usage if (@ARGV !=4 ||!$chr_id);

$r_type ||='first';
$name ||='test';
$outdir ||=getcwd();

print 'Start time: '.localtime()."\n";
my @solution_set;
if (-e $path){
    open INA,$path or die $!;
    while(<INA>){
        chomp;
        my $solu = (split/\s+/,$_)[1];
        push @solution_set,$solu;
    }
}else{
    push @solution_set,$path;
}

my $min_edge;
open INB,$edge_file or die $!;
while(<INB>){
    chomp;
    next if (/^\#/);
    my @set = split/\s+/,$_;
    next if ($set[1] == 0);
    $min_edge = $set[0];
    last;
    
}
close INB;

my %edge_region_set = mersh_edge_region($edge_file,$region_file);
open OUTB,">$outdir/$name.path.index" or die $!;
if ($r_type eq 'all'){
    my $num = 1;
    foreach my $solution (@solution_set){
        open OUTA,">$outdir/$name.$num.fa" or die $!;
        my $new_path = $solution;
        if ($chr_type eq 'cycle'){
            my @edge_set = split/\//,$solution;
            my @min_edge_index = indexes{$_ eq $min_edge} @edge_set;
            unless (@min_edge_index){
                my $rever_path = link_reverse($solution);
                @edge_set = split/\//,$rever_path;
                @min_edge_index = indexes {$_ eq $min_edge} @edge_set;
            }
            my $element = scalar(@edge_set);
            my @first_set = @edge_set[0..$min_edge_index[-1]-1];
            my @second_set = @edge_set[$min_edge_index[-1]..$element-1];
            my @new_set = (@second_set,@first_set);
            $new_path = join("/",@new_set);
        }
        my @path_set = path_index($new_path,\%edge_region_set);
        print OUTB "$name\t$new_path\t".join(",",@path_set)."\n";
        my $this_fasta = translate_path($fasta,$new_path,$chr_id);
        my $length = length($this_fasta);
        print OUTA ">$chr_id#$num#$length"."\n";
        for (my $i=0;$i<$length;$i=$i+60){
            my $line_fa = substr($this_fasta,$i,60);
            print OUTA $line_fa."\n";
        }
        close OUTA;
        $num++;
    }
}elsif($r_type eq 'first'){
    open OUTA,">$outdir/$name.1.fa" or die $!;
    #my @index_set = path_index($path_set[0],\%edge_region_set);
    my $new_path = $solution_set[0];
    if ($chr_type eq 'cycle'){
        my @edge_set = split/\//,$solution_set[0];
        my @min_edge_index = indexes{$_ eq $min_edge} @edge_set;
        unless (@min_edge_index){
            my $rever_path = link_reverse($solution_set[0]);
            @edge_set = split/\//,$rever_path;
            @min_edge_index = indexes {$_ eq $min_edge} @edge_set;
        }
        my $element = scalar(@edge_set);
        my @first_set = @edge_set[0..$min_edge_index[-1]-1];
        my @second_set = @edge_set[$min_edge_index[-1]..$element-1];
        my @new_set = (@second_set,@first_set);
        $new_path = join("/",@new_set);
    }
    my @path_set = path_index($new_path,\%edge_region_set);
    print OUTB "$name\t$new_path\t".join(",",@path_set)."\n";
    my $this_fasta = translate_path($fasta,$new_path,$chr_id);
    my $length = length($this_fasta);
    print OUTA ">$chr_id#1#$length"."\n";
    for (my $i=0;$i<$length;$i=$i+60){
        my $line_fa = substr($this_fasta,$i,60);
        print OUTA $line_fa."\n";
    }
    close OUTA;
}else{
    print 'please input the corret paramaters !!!'."\n";
}
close OUTB;
print 'End time: '.localtime()."\n";
print '######## This work is finished! #########'."\n";
#################################################################
sub translate_path{
    my ($fasta,$path,$chrid) = @_;
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
            $rev_seq =~ tr/ATCGatcg/TAGCtagc/;
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

###########################################
sub path_index{
    my ($input,$edge_hash)=@_;
    my @list = split/\//,$input;
    my @output_set;
    foreach my $aedge (@list){
        my ($lindex,$rindex) = (split/_/,$aedge)[0,1];
        if ($lindex > $rindex){
            my $bedge = $rindex.'_'.$lindex;
            my $index1 = ${${$edge_hash}{$bedge}}[4];
            my @index_set = split/\,/,$index1;
            my @n_set;
            foreach my $n (@index_set){
                my $new_n = '-'.$n;
                unshift @n_set,$new_n;  
            }
            push @output_set,@n_set;
        }else{
            my $index2 = ${${$edge_hash}{$aedge}}[4];
            push @output_set,$index2;
        }
        
    }
    return @output_set;
}
###################################
sub link_reverse{
    my $edge2 = $_[0];
    my @set1 = split/\//,$edge2;
    my @new_set;
    foreach my $ed (@set1){
        my ($l,$r) = split/_/,$ed;
        my $new_ed = $r."_".$l;
        unshift @new_set,$new_ed;
    }
    my $new_result = join("/",@new_set);
    return $new_result;
}

