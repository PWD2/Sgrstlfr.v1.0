#!/usr/bin/perl -w

use strict;
<<<<<<< HEAD
#use Data::Printer;
=======
use Data::Printer;
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
use Cwd qw/getcwd/;
use List::Util qw/max/;
use List::MoreUtils qw/uniq/;
use Array::Utils qw(:all);

my $usage=<<USAGE;

    Descripton: v1.0 2023.01.04 pwd
    This script is order to class the node barcode by the basic edge and nod inforamtion!

    Example: perl $0 <nodebarcode> <notnodebarcode> <edge> <nod> <name> <outdir>

            <nodebarcode>    : the barcode information have node
            <notnodebarcode> : the barcode information not have node
            <edge>           : the edge inforamtion
            <nod>            : the node information
            <name>           : output file prefix
            <outdir>         : output directory

USAGE

my ($b1,$b2,$edge,$nod,$name,$outdir) = @ARGV;

die $usage if (!$b1||!$b2||!$edge||!$nod);

$name ||='test';
$outdir ||=getcwd;

print 'Start time '.localtime()."\n";
my %edge_cn_hash;
my @mutiple_cn_set;
open INB,$edge or die $!;
while(<INB>){
    chomp;
    next if (/^\#/);
    my ($edg,$cn) = (split,/\s+/,$_)[0,1];
    next if ($cn == 0);
    if ($cn == 1){
        $edge_cn_hash{$edg} = $cn;
    }else{
        push @mutiple_cn_set,$edg;
    }
    
}
close INB;

my %nod_all;
open INC,$nod or die $!;
while(<INC>){
    my ($a_n,$num) = split/\s+/,$_;
    my ($le,$ri) = split/\//,$a_n;
    push @{$nod_all{$le}},($a_n);
    my $ri_an = $ri.'/'.$le;
    push @{$nod_all{$ri}},($ri_an);
}
close(INC);

### repiar the nodebarcode set
$/="\@";
my %sort_barcode_set;
my %merge_barcode_set;
my $node_all_num=0;
my $node_have_uniq_edge_num=0;
my $node_nothave_uniq_edge_num=0;
my $node_overlap_num=0;
my $sort_barcode_stat=0;
open OUTA,">$outdir/$name.rn.barcode" or die $!;
open OUTB,">$outdir/$name.merge.barcode" or die $!;
open OUTD,"| sort -k 2 -k 4 -n >$outdir/$name.sort.barcode.stat" or die $!;
open INA,$b1 or die $!;
<INA>;
while(<INA>){
    chomp;
    my ($id,$dis,$bili,$node) = split/\n/,$_;
    my ($barcode,$edge_nod) = split/\s+/,$id;
    my @edge_set = split/\,/,$edge_nod;
    my $this_edge_num = scalar(@edge_set);
    my @node_set = split/\,/,$node;
    my $this_node_num = scalar(@node_set);
    $node_all_num++;
    ####repair the nod set
    my @this_uniq_edge;
    foreach my $aedge (@edge_set){
        if (exists $edge_cn_hash{$aedge}){
            push @this_uniq_edge,$aedge;
        }
    }
    if (@this_uniq_edge){
        my @new_nod_set;
        my @include_nod;
        $node_have_uniq_edge_num++;
        foreach my $u_edge (@this_uniq_edge){
            my ($l,$r)=split/_/,$u_edge;
            if (exists $nod_all{$l}){
                push @new_nod_set,${$nod_all{$l}}[0];
            }
            if (exists $nod_all{$r}){
                push @new_nod_set,${$nod_all{$r}}[0];
            }
        } 
        my @uniq_new_nod_set = uniq(re_sort_nod(\@new_nod_set));
        push @include_nod,@uniq_new_nod_set;
        ###repair the edge set
        my @repair_edge_set = ntoedge(\@node_set,'a');
        my @new_edge_set = uniq(@edge_set,@repair_edge_set);
        my $new_edge_scalar = scalar(@new_edge_set);
        ###merge the result
        my @the_nod_set = uniq(@include_nod,@node_set);
        my %nod_edge = ntoedge(\@the_nod_set,'h');
        my $over_num = 0 ;
        foreach my $this_key (keys %nod_edge){
            my @this_nod_set = @{$nod_edge{$this_key}};
            if (intersect(@edge_set,@this_nod_set) == 2){
                $over_num++;
            }
            else{
               delete $nod_edge{$this_key};
            }
        }
        my @final_nod_set = keys %nod_edge;
        my $re_edge = join(",",@new_edge_set);
        my $re_edge_num = scalar(@new_edge_set);
        my $re_node = join(",",@final_nod_set);
        my $re_nod_num = scalar(@final_nod_set);
        print OUTB "$barcode\t$re_edge_num\t$re_edge\t$re_nod_num\t$re_node\t$bili\n";
        if ($over_num >= $new_edge_scalar - 1){
            print OUTA '@'."$barcode\t".join(",",@new_edge_set)."\n$dis\n$bili\n".join(",",@final_nod_set)."\n";
            $node_overlap_num++;
            my $edge_number = scalar(@new_edge_set);
            my $nod_number = scalar(@final_nod_set);
            if ($edge_number - $nod_number <=1){
                print OUTD "$barcode\t$edge_number\t".join(",",@new_edge_set)."\t$nod_number\t".join(",",@final_nod_set)."\t$bili\n";
                $sort_barcode_stat++;
            }
        }
    }
    else{
        $node_nothave_uniq_edge_num++;
    }
}
close INA;

print 'The node barcode set are all of number: '.$node_all_num."\n";
print 'The node barcode set have uniq edge number: '.$node_have_uniq_edge_num."\n";
print 'The node barcode set have all duplication edge number: '.$node_nothave_uniq_edge_num."\n";
print 'The node barcode set have overlap number: '.$node_overlap_num."\n\n";
##repair the notnod barcode set 
my $notnode_all_num=0;
my $notnode_dup_edge_num=0;
my $notnode_uniq_edge_num=0;
my $notnode_overlap_num=0;
open OUTC,">$outdir/$name.nd.barcode" or die $!;
open IND,$b2 or die $!;
<IND>;
while(<IND>){
    chomp;
    my ($id,$dis,$bili,$node) = split/\n/,$_;
    my ($barcode,$edge_nod) = split/\s+/,$id;
    my @edge_set = split/\,/,$edge_nod;
    my $this_edge_num = scalar(@edge_set);
    my $number=0;
    my @this_uniq_edge;
    my @mutiple_cn_edge;
    $notnode_all_num++;
    foreach my $aedge (@edge_set){
        if (exists $edge_cn_hash{$aedge}){
            $number++;
            push @this_uniq_edge,$aedge;
        }else{
            push @mutiple_cn_edge,$aedge;
        }
    }
    if ($number == 0){
        print OUTC $_;
        $notnode_dup_edge_num++;
    }else{
        $notnode_uniq_edge_num++;
        my @new_nod_set;
        foreach my $u_edge (@this_uniq_edge){
            my ($l,$r)=split/_/,$u_edge;
            if (exists $nod_all{$l}){
                push @new_nod_set,${$nod_all{$l}}[0];
            }
            if (exists $nod_all{$r}){
                push @new_nod_set,${$nod_all{$r}}[0];
            }
        } 
        my @uniq_new_nod_set = uniq(re_sort_nod(\@new_nod_set));
        my %nod_edge = ntoedge(\@uniq_new_nod_set,'h');
        my $over_num = 0 ;
        foreach my $this_key (keys %nod_edge){
            my @this_nod_set = @{$nod_edge{$this_key}};
            if (intersect(@edge_set,@this_nod_set) == 2){
                $over_num++;
            }
            else{
               delete $nod_edge{$this_key};
            }
        }
        my @final_nod_set = keys %nod_edge;
        if (@final_nod_set == 0){
            print OUTC $_;
            $notnode_dup_edge_num++;
        }else{
            my $re_node = join(",",@final_nod_set);
            my $re_node_num = scalar(@final_nod_set);
            print OUTB "$barcode\t$this_edge_num\t$edge_nod\t$re_node_num\t$re_node\t$bili\n";
            if ($over_num >= $this_edge_num - 1){
                print OUTA '@'."$id\n$dis\n$bili\n".join(",",@final_nod_set)."\n";
                $notnode_overlap_num++;
                my ($bcode,$t_edge) = split/\s+/,$id;
                my $t_edge_number = split/\,/,$t_edge;
                my $t_nod_number = scalar(@final_nod_set);
                if ($t_edge_number - $t_nod_number <= 1){
                    print OUTD "$bcode\t$t_edge_number\t$t_edge\t$t_nod_number\t".join(",",@final_nod_set)."\t$bili\n";
                    $sort_barcode_stat++;
                }
            }
        }
        
    }
}
$/="\n";
close IND;

close OUTC;
close OUTD;

print 'The not node barcode set are all of number: '.$notnode_all_num."\n";
print 'The not node barcode set have all duplication edge number: '.$notnode_dup_edge_num."\n";
print 'The not node barcode set have uniq edge number: '.$notnode_uniq_edge_num."\n";
print 'The not node barcode haver overlap number: '.$notnode_overlap_num."\n\n";

print 'The all barcode number: '.($notnode_all_num+$node_all_num)."\n";
print 'We get the barcode have uniq edge number: '.($notnode_uniq_edge_num+$node_have_uniq_edge_num)."\n";
print 'We get the barcode have overlap number: '.($notnode_overlap_num+$node_overlap_num)."\n";
print 'We get the sort barcode set number: '.$sort_barcode_stat."\n";

print 'End time: '.localtime()."\n";
######################################################
sub re_sort_nod{
    my $input =$_[0];
    my @this_set;
    if (ref($input) eq 'ARRAY'){
        @this_set = @$input;
    }else{
        @this_set = split/_/,$input;
        pop @this_set;
        shift @this_set;
    }
    my @new_set;
    foreach my $this_nod (@this_set){
        my ($l,$r) = split/\//,$this_nod;
        if ($l > $r){
            my $a_nod = $r.'/'.$l;
            push @new_set,$a_nod;
        }else{
            push @new_set,$this_nod;
        }
    }
    return @new_set;
}

#####################################
sub ntoedge{
    my ($input,$type) = @_;
    if ($type eq 'h'){
        my %a_nod_set;
        foreach my $tnod (@$input){
            my ($l,$r) = split/\//,$tnod;
            my $left_sight = judge_odd_even($l,'l');
            my $right_sight = judge_odd_even($r,'r');
            my @this_set = ($left_sight,$right_sight);
            my @sort_edge = re_sort_edge(\@this_set);
            @{$a_nod_set{$tnod}} = @sort_edge;
        }
        return %a_nod_set;
    }else{
        my @a_nod_set;
        foreach my $tnod (@$input){
            my ($l,$r) = split/\//,$tnod;
            my $left_sight = judge_odd_even($l,'l');
            my $right_sight = judge_odd_even($r,'r');
            push @a_nod_set,($left_sight,$right_sight);
        }
        my @sort_nod_set = re_sort_edge(\@a_nod_set);
        return @sort_nod_set;
    }
}

############################################
sub re_sort_edge{
    my $input = $_[0];
    my @new_set;
    foreach my $edge (@$input){
        my ($l,$r)=(split/_/,$edge)[0,1];
        if ($l > $r){
            my $new_edge = $r.'_'.$l;
            push @new_set,$new_edge;
        }else{
            push @new_set,$edge;
        }
    }
    return @new_set;
}

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

sub rm_mutiple_nod_solution{
    my $input = $_[0];
    my %number_stat;
    foreach my $apn (@$input){
        my ($l,$r) = split/\//,$apn;
        push @{$number_stat{$l}},($apn);
        push @{$number_stat{$r}},($apn);
    }
    my $st_num=0;
    foreach my $tky (keys %number_stat){
        my @set = @{$number_stat{$tky}};
        my $s_n = scalar(@set);
        if ($s_n >1){
            $st_num++;
        }
    }
    if ($st_num > 1){
        return 'F';
    }else{
        return 'T';
    }
}

