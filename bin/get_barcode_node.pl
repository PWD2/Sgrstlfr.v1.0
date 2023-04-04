#!/usr/bin/perl -w

use strict;
use Data::Printer;
use Basicread qw/mersh_edge_region/;
use List::Util qw/pairs/;
use List::MoreUtils qw/uniq/;
use POSIX qw/floor/;
use Getopt::Long;
use Cwd qw/getcwd/;

my ($edge,$nod,$loxpregion,$map,$synmap) = @ARGV;
my ($chrtype,$name,$insert_size,$overlap,$outdir);
my $usage=<<USAGE;

Description: v1.0 2023.02.13 pwd 
             This script is purpose to get the barcode nod by the PE and SR and get
             the barcode other complete information!!!

Usage: perl $0 <edge> <nod> <loxpregion> <split.sam> <synmap> [-option]

            -t              :[linear | cycle] the chrmosome type [must] 
            -isz            :the insert size of the pair-end seqence [default:500]
            -overlap        :the overlap of reads in range! [default:0.8]
            -n              :the file name
            -o              :the output directory

USAGE

GetOptions(
    "t:s" => \$chrtype,
    "isz:s" => \$insert_size,
    "overlap:s" =>\$overlap,
    "n:s" => \$name,
    "o:s" => \$outdir
);

die $usage if (!$edge||!$nod||!$map||!$synmap||!$loxpregion||!$chrtype);

$chrtype||='linear';
$name ||='test';
$overlap ||= 0.8;
$insert_size ||= 500;
$outdir ||=getcwd;
my $loxpsym = 'ATAACTTCGTATAATGTACATTATACGAAGTTAT';

print '#Start:'.localtime()."\n";
print "A1: Start to get the barcode node, loading.......\n";
###########Step1 make sure the range of the number ##################
my %edge_region = mersh_edge_region($edge,$loxpregion);
my %repair_hash;
my %range_hash;
foreach my $keys (keys %edge_region){
    my @set = @{$edge_region{$keys}};
    next if ($set[0] == 0);
    my ($left_ldex,$right_ldex) = split/_/,$keys;
    my $median = $set[1]/2;
    my $plus1 = 110;
    if ($median < 110){
        $plus1 = $median;
    }
    @{$repair_hash{$left_ldex}} = ($set[2],$set[2]+$plus1 - 1);
    @{$repair_hash{$right_ldex}} = ($set[3]-$plus1 + 1,$set[3]);  
    
    if ($set[1] < $insert_size){
        next;
    }else{
        my $plus2 = $insert_size;
        if ($insert_size > $median){
            $plus2 = $median;
        }
        push @{$range_hash{$left_ldex}},($set[2],$set[2] + $plus2 -1) ;   
        push @{$range_hash{$right_ldex}},($set[3] - $plus2 + 1,$set[3]) ;    
    }
}

#########Step2 get the true start_point by loxpsym ###########################
open INB,$map or die $!;
my %pair_reads;
while(<INB>){
    chomp;
    next if (/^\@/);
    my ($id,$start_point1,$my_seq) = (split/\s+/,$_)[0,3,9];
    next if ($start_point1 =~/^\*/);
    my $check_start;
    if ($my_seq =~ /^$loxpsym/gi){
        $check_start = $start_point1 + 17;
    }elsif($my_seq =~ /$loxpsym$/gi){
        $check_start = $start_point1;
    }else{
        print "$id\n";
        print '!why not match the loxpsym???'."\n";
    }
    push @{$pair_reads{$id}},($check_start);
}
close INB;
# p(%pair_reads);
# exit;

###############Step3 locate the reads place in the range nod #################
#SR-nod
my %barcode_nod;
foreach my $keys1 (keys %pair_reads){
    my @reads_set = uniq(@{$pair_reads{$keys1}});
    my $pair_num = scalar(@reads_set);
    next if ($pair_num % 2 == 1);
    my $barcode_id;
    if ($keys1 =~ /\#(.+)/){
        $barcode_id = $1;
    }
    if ($pair_num <= 3){
        my $a_l = corret_pe($reads_set[0],\%repair_hash);
        my $a_r = corret_pe($reads_set[1],\%repair_hash);
        my $this_nod = make_sure_pe($a_l,$a_r,$chrtype);
        push @{$barcode_nod{$barcode_id}},($this_nod) if ($this_nod ne 'undef');
    }else{
        foreach my $pair (pairs @reads_set){ 
            my ($key,$value) = @$pair;
            my $a_l = corret_pe($key,\%repair_hash);
            my $a_r = corret_pe($value,\%repair_hash);
            my $this_nod = make_sure_pe($a_l,$a_r,$chrtype);
            push @{$barcode_nod{$barcode_id}},($this_nod) if ($this_nod ne 'undef');        
        }
    }
}

#PE-nod
my %syn_map;
my $readsnum=0;
if ($synmap =~ /\.gz$/){
    open INC,"gzip -dc $synmap |" or die $!;
}else{
    open INC,"$synmap" or die $!;
}
while(<INC>){
    chomp;
    my ($id,$start_point1,$end_point1,$start_point2,$end_point2) = (split/\s+/,$_)[0,2,3,6,7];
    my $barcode;
    if ($id =~ /\#(.+)/){
        $barcode = $1;
    }
    if ($barcode ne '0_0_0'){
        push @{$syn_map{$barcode}},[$start_point1,$end_point1];
        push @{$syn_map{$barcode}},[$start_point2,$end_point2];
        $readsnum = $readsnum + 2;
        my $id_l = pe_belong_range($start_point1,$end_point1,$overlap,\%range_hash);
        my $id_r = pe_belong_range($start_point2,$end_point2,$overlap,\%range_hash);
        my $t_nod = make_sure_pe($id_l,$id_r,$chrtype);
        push @{$barcode_nod{$barcode}},($t_nod) if ($t_nod ne 'undef');
    }
}
close INC;
print '#The useful reads number :'."$readsnum\n";
##########Step4 check the barcode nod if corret ############################
my %check_nod;
my @corret_set;
open IND,$nod or die $!;
while(<IND>){
    chomp;
    my $node = (split/\s+/,$_)[0];
    push @corret_set,$node;
}
close IND;

my $all_nod_barcode=0;
my $one_nod_barcode=0;
my $two_nod_barcode=0;
my $more_than_two=0;

open OUT,">$outdir/$name.barcode.nod" or die $!;
foreach my $keys2 (keys %barcode_nod){
    my @this_nod = @{$barcode_nod{$keys2}};
    my @a_new_set;
    foreach my $keys1 (@this_nod){
        my ($left,$right) = split/\//,$keys1;
        my $new_node;
        if ($left > $right){
            $new_node = $right.'/'.$left;
        }else{
            $new_node = $keys1
        }
        if (grep {$new_node eq $_} @corret_set){
            push @a_new_set,$new_node;
        }
        # else the erro nod
    }
    my @uniq_set = uniq(@a_new_set);
    if($uniq_set[0]){
        $all_nod_barcode++;
        if (@uniq_set == 1){
            $one_nod_barcode++;
        }
        if (@uniq_set == 2){
            $two_nod_barcode++;
        }
        if (@uniq_set > 2){
            $more_than_two++;
        }
        print OUT "$keys2\t".join(",",@uniq_set)."\n";
        push @{$check_nod{$keys2}},(@uniq_set);
    }
}
close OUT;

print '#The barcode have nod number: '.$all_nod_barcode."\n";
print '#The one nod barcode: '.$one_nod_barcode."\n";
print '#The two nod barcode: '.$two_nod_barcode."\n";
print '#The more than two nod barcode: '.$more_than_two."\n";

print 'A1 is done !'."\n";

print 'A2 Start to get the barcode edge !!'."\n";
#########Step5 get the place of reads in edge range ##################
my $statreads=0;
my %stat_hash;
open OUTC,">$outdir/$name.co_barcode_reads_freq.txt" or die $!; 
foreach my $barcode (keys %syn_map){
    my @reads_set = @{$syn_map{$barcode}};
    my $reads_num = scalar(@reads_set);
    print OUTC "$barcode\t$reads_num\n";
    foreach my $reads (@reads_set){
        my $read1_start = ${$reads}[0];
        my $read1_end = ${$reads}[1];
        my ($p1_edge_type,$r_start) = reads_belong_edge($read1_start,$read1_end,$overlap,\%edge_region);       
        if ($p1_edge_type ne 'undef'){
            push @{$stat_hash{$barcode}{$p1_edge_type}},($r_start);
            $statreads = $statreads + 1; 
        }else{
            next;
        } 
    }
   
}
print "#The useful reads number: $statreads\n";
close OUTC;
#########Step6 stat the rpkm and collect the all inforamtion #############
open OUTA,">$outdir/$name.node.barcode" or die $!;
open OUTB,">$outdir/$name.notnode.barcode" or die $!;

my $max_edge_number = 0;
foreach my $this_barcode (keys %stat_hash){
    my %second_hash = %{$stat_hash{$this_barcode}};
    #my $reads_num = keys %second_hash;
    next if (keys %second_hash < 3);
    my %rpkm_edge;
    my @sort_set;
    my @print_set;
    my $all_reads_num = 0;
    foreach my $this_edge ( sort {(split/_/,$a)[0]<=>(split/_/,$b)[0]} keys %second_hash){
        my @read_set = @{$second_hash{$this_edge}};
        my $this_number = scalar(@read_set);
        $all_reads_num = $all_reads_num + $this_number;
        my $this_length = ${$edge_region{$this_edge}}[1];
        my $this_cn = ${$edge_region{$this_edge}}[0];
        ###rpkm ####
        push @{$rpkm_edge{$this_edge}},($this_cn,$this_number,$this_length);
        ###mapping point line ####
        my $output_line = $this_edge.';'.join(",",sort {$a<=>$b} @read_set).';'.$this_number;
        push @print_set,$output_line;
        push @sort_set,$this_edge;
    }
    my $average_edge_reads = $all_reads_num/keys %second_hash;
    #next if ($average_edge_reads < $filter);
    my $edge_number = scalar(@sort_set);
    if ($edge_number > $max_edge_number){
        $max_edge_number = $edge_number;
    }
    if ($check_nod{$this_barcode}){
        ###line 1 ####    
        print OUTA '@'.$this_barcode."\t".join(",",@sort_set)."\n";
        ###line 2 ####
        print OUTA join("\t",@print_set)."\n";
        ###line 3 ####
        my %rpkm_set = stat_rpkm(\%rpkm_edge,$statreads);
        my $k=1;
        foreach my $keys4 (sort {(split/_/,$a)[0]<=>(split/_/,$b)[0]} keys %rpkm_set){
            my $bi_li = $rpkm_set{$keys4};
            print OUTA ';' if ($k > 1);
            print OUTA $keys4.':'.$bi_li;
            $k++;
        }
        print OUTA "\n";
        ###line 4 ####
        print OUTA join(",",@{$check_nod{$this_barcode}})."\n";
    }else{
        ###line 1 ####    
        print OUTB '@'.$this_barcode."\t".join(",",@sort_set)."\n";
        ###line 2 ####
        print OUTB join("\t",@print_set)."\n";
        ###line 3 ####
        my %rpkm_set = stat_rpkm(\%rpkm_edge,$statreads);
        my $k=1;
        foreach my $keys4 (sort {(split/_/,$a)[0]<=>(split/_/,$b)[0]} keys %rpkm_set){
            my $bi_li = $rpkm_set{$keys4};
            print OUTB ';' if ($k > 1);
            print OUTB $keys4.':'.$bi_li;
            $k++;
        }
        print OUTB "\n";
        ###line 4 ####
        print OUTB '+'."\n";
    }
}

close OUTA;
close OUTB;

print 'A2 is done !'."\n";
print '#The max edge number is '.$max_edge_number."\n";
print 'This work is finshed!'."\n";
print '#End:'.localtime()."\n";
##########################################
sub make_sure_pe{
    my ($a,$b,$type) = @_;
    my $a_nod = 'undef';
    if ($a ne 'undef' and $b ne 'undef'){
        if ($type eq 'liner'){
            if ($a != 0 and $b !=0){
                if ($a != $b){
                    if ($a < $b){
                        $a_nod = $a.'/'.$b;
                    }else{
                        $a_nod = $b.'/'.$a;
                    }
                }
            }
        }else{
            if ($a != $b){
                if ($a<$b){
                    $a_nod = $a.'/'.$b;
                }else{
                    $a_nod = $b.'/'.$a;
                }
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
    my $result_line = 'undef';
    my $repair_start ='undef';
    foreach my $line (keys %$range_hash){
        my $range_a = ${${$range_hash}{$line}}[2];
        my $range_b = ${${$range_hash}{$line}}[3];
        ##contain
        if ($range_a <= $read_a and $read_b <= $range_b){
            $result_line = $line;
            $repair_start = $read_a - $range_a;
        ##left sight
        }elsif($read_a < $range_a and $read_b > $range_a){
            my $precent = ($read_b - $range_a + 1)/100;
            if ($precent >= $include){
                $result_line = $line;
                $repair_start = $read_a - $range_a;
            }
        ##right sight
        }elsif($read_a < $range_b and $range_b < $read_b){
            my $precent = (100-($read_b - $range_b + 1))/100;
            if ($precent >= $include){
                $result_line = $line;
                $repair_start = $read_a - $range_a;
            }
        }else{
            next;
        }
    }
    return ($result_line,$repair_start);
}
