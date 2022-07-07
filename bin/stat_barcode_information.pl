#!/usr/bin/perl -w 

use strict;
use Data::Printer;
use Cwd;
use Basicread qw/mersh_edge_region/;
use List::Util qw/max min/;
use List::MoreUtils qw/uniq/;
use Statistics::Descriptive;
use Getopt::Long;

my $usage=<<USAGE;

Dsecribe: This script is purpose to get the complete barcode information of co-barcode reads.

          Author:pangwending\@genomics.cn;
          Version: 1.0, Date: 2022.06.09

Usage:perl $0 <edge> <loxpregion> <synmap> <barcode_node>[-option] (edge file format:line,num,depeth,region)

            -n          : Sample name [default:test]
            -f          : filter the barcode by the per-edge average reads number [defacult:4]
            -overlap    : the overlap of the reads in edge_region(0-1) [default:0.8]
            -o          : the output directory [defacult:./]
            -h          : show this help
            
Example: perl $0 JS599.edg IXR_BACseq.loxpreg JS599.synmap.gz [-option].

USAGE
my ($edge_list,$loxplist,$synmap,$bsnode)=@ARGV;
my ($name,$overlap,$filter,$outdir,$help);
GetOptions(
    "n:s"=> \$name,
    "f:s" => \$filter,
    "overlap:s"=> \$overlap,
    "o:s"=> \$outdir,
    "h:s"=> \$help
);

die $usage if (!$edge_list||!$loxplist||!$synmap||$help);

$name ||='test';
$filter ||='4';
$overlap ||='0.8';
$outdir ||=getcwd;

print '#Start:'.localtime()."\n";
print "#output the sub-edge, loading.......\n";

############## get the edge and region set by the modelu ##################################
my %edge_region = mersh_edge_region($edge_list,$loxplist);
foreach my $keys4 (keys %edge_region){
    my @a_set = @{$edge_region{$keys4}};
    delete $edge_region{$keys4} if ($a_set[0] == 0);
}

###############Step1 get the pair-reads information #######################################
if ($synmap =~ /\.gz$/){
    open INC,"gzip -dc $synmap |" or die $!;
}else{
    open INC,$synmap || die $!;
}
my %syn_map;
my $readsnum=0;
while(<INC>){
    chomp;
    my ($id,$start1,$end1,$start2,$end2)=(split/\s+/,$_)[0,2,3,6,7];  #pair reads1 and reads 2 
    if ($id =~ /\#(.+)/){
        my $barcode = $1;
        push @{$syn_map{$barcode}},[$start1,$end1] if($barcode ne '0_0_0');
        push @{$syn_map{$barcode}},[$start2,$end2] if($barcode ne '0_0_0');
        $readsnum = $readsnum + 2;
    }else{
        print "This $id reads no barcode???\n";
    }
}
#p(%syn_map);exit;
close INC;
print "#reads number:$readsnum\n";

##############Step2 get the mersh node set(include split-nod and pe-nod) ########################
open IND,$bsnode or die $!;
my %barcode_node_set;
while(<IND>){
    chomp;
    my ($bs,$nod) = split/\s+/,$_;
    $barcode_node_set{$bs} = $nod;
}
close IND;

#############Step3 get the place of reads in the edge #####################################
my $statreads=0;
my %stat_hash; 
my $max_edge_num =1;
foreach my $barcode (keys %syn_map){
    my @reads_set = @{$syn_map{$barcode}};
    foreach my $reads (@reads_set){
        my $read1_start = ${$reads}[0];
        my $read1_end = ${$reads}[1];
        my $p1_edge_type = reads_belong_edge($read1_start,$read1_end,$overlap,\%edge_region);       
        if ($p1_edge_type ne 'undef'){
            push @{$stat_hash{$barcode}{$p1_edge_type}},($read1_start);
            $statreads = $statreads + 1; 
        }else{
            next;
        } 
    }
   
}
print "#The useful reads number: $statreads\n\n";

############Step4 stat the tmp to make sure the prpportion and collect and all result ###############
open OUTA,">$outdir/$name.node.barcode" or die $!;
open OUTB,">$outdir/$name.notnode.barcode" or die $!;
foreach my $keys1 (keys %stat_hash){
    my %second_hash = %{$stat_hash{$keys1}};
    next if (keys %second_hash < 3);
    my %rpkm_edge;
    my @sort_set;
    my @print_set;
    my $all_reads_num = 0;
    foreach my $keys2 ( sort {(split/_/,$a)[0]<=>(split/_/,$b)[0]} keys %second_hash){
        my @read_set = @{$second_hash{$keys2}};
        my $this_number = scalar(@read_set);
        $all_reads_num = $all_reads_num + $this_number;
        my $this_length = ${$edge_region{$keys2}}[1];
        my $this_cn = ${$edge_region{$keys2}}[0];
        ###rpkm ####
        push @{$rpkm_edge{$keys2}},($this_cn,$this_number,$this_length);
        ###mapping point line ####
        my $output_line = $keys2.';'.join(",",sort {$a<=>$b} @read_set).';'.$this_number;
        push @print_set,$output_line;
        push @sort_set,$keys2;
    }
    my $average_edge_reads = $all_reads_num/keys %second_hash;
    next if ($average_edge_reads < $filter);
    if ($barcode_node_set{$keys1}){
        ###line 1 ####    
        print OUTA '@'.$keys1."\t".join(",",@sort_set)."\n";
        ###line 2 ####
        print OUTA join("\t",@print_set)."\n";
        ###line 3 ####
        my %rpkm_set = stat_rpkm(\%rpkm_edge,$statreads);
        my $k=1;
        foreach my $keys3 (sort {(split/_/,$a)[0]<=>(split/_/,$b)[0]} keys %rpkm_set){
            my $bi_li = $rpkm_set{$keys3};
            print OUTA ';' if ($k > 1);
            print OUTA $keys3.':'.$bi_li;
            $k++;
        }
        print OUTA "\n";
        ###line 4 ####
        print OUTA $barcode_node_set{$keys1}."\n";
    }else{
        ###line 1 ####    
        print OUTB '@'.$keys1."\t".join(",",@sort_set)."\n";
        ###line 2 ####
        print OUTB join("\t",@print_set)."\n";
        ###line 3 ####
        my %rpkm_set = stat_rpkm(\%rpkm_edge,$statreads);
        my $k=1;
        foreach my $keys3 (sort {(split/_/,$a)[0]<=>(split/_/,$b)[0]} keys %rpkm_set){
            my $bi_li = $rpkm_set{$keys3};
            print OUTB ';' if ($k > 1);
            print OUTB $keys3.':'.$bi_li;
            $k++;
        }
        print OUTB "\n";
        ###line 4 ####
        print OUTB '+'."\n";
    }
}

close OUTA;
close OUTB;
#################stat rpkm by complete genome ###################
sub stat_rpkm{
    my ($set,$total_reads) = @_;
    my %new_hash;
    my $total_rpkm = 0;
    foreach my $keys (keys %$set){
        my @this_set = @{${$set}{$keys}};
        my $first_stat = $this_set[1]*10000000/$total_reads;
        my $rpkm = $first_stat/$this_set[2];
        $total_rpkm = $total_rpkm + $rpkm;
        $new_hash{$keys} = $rpkm;
    }
    my $edge_num = keys %new_hash;
    my $average = $total_rpkm/$edge_num;
    my %bi_li_edge;
    foreach my $keys1 (keys %new_hash){
        my $a_rpkm = $new_hash{$keys1};
        my $acn = ${${$set}{$keys1}}[0];
        for (my $i=0;$i<$acn;$i++){
            if ($i*$average < $a_rpkm and $a_rpkm <= ($i+1)*$average){
                $bi_li_edge{$keys1} = ($i+1);
            }else{
                next;
            }
        }
        unless(exists $bi_li_edge{$keys1} ){
            $bi_li_edge{$keys1} = $acn;
        }
    }
    return %bi_li_edge;
}

##################sort_edge ######################
sub sort_edge{
    my $input = $_[0];
    my %hash;
    foreach my $edge (@$input){
        my $l = (split/_/,$edge)[0];
        $hash{$l} = $edge;
    }
    my @new;
    foreach my $keys (sort {$a <=> $b} keys %hash){
        push @new,$hash{$keys};
    }
    return @new;
}

####################################################
sub reads_belong_edge{
    my ($read_a,$read_b,$include,$range_hash) = @_;
    my $result_line;
    foreach my $line (keys %$range_hash){
        my $range_a = ${${$range_hash}{$line}}[2];
        my $range_b = ${${$range_hash}{$line}}[3];
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


