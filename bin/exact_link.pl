#!/usr/bin/perl -w

use strict;
<<<<<<< HEAD
#use Data::Printer;
=======
use Data::Printer;
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
use Cwd qw/getcwd/;
use Getopt::Long;
use threads;
use List::Util qw/max min sum tail head/;
use Array::Utils qw(:all);
use List::MoreUtils qw/uniq/;
use Classbarcode qw/class_barcode same_class_barcode/;
use POSIX qw/floor/;

my $usage=<<USAGE;

    Description: This script is order to aseembly the SCRaMbLE genome by the LFR which from stLFR seqence !!!!

    v1.0 pangwending 2023.01.04
    perl $0 <sort.barcode.stat> <nod> <edg> <total> [-option]

        -n              :Sample name
        -o              :outdir
<<<<<<< HEAD
        -chr_type       :[linear | cycle] the chrmosome type
=======
        -chr_type       :[liner | cycle] the chrmosome type
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
        -cf             :the genome erro(confict) about the edge and nod [default:0] 
        -fix            :[undef | others]Please input the new duplication type! [example: 14_15,T,T;8_9,T,F (edge,tandem,untandem)]
        -sf             :the fix scaffold set as input [defult:undef]
        -b              :the link-edge barcode use max [N+b] [defult:b=1]
        -fp             :[auto | others]the number of the link-path barcode support. [defult:auto]
        -t              :the threads of more than 5000 solution will use [defult:4] 
        -i              :the allow erro in the scaffold set [default:0]
        -len            :the start cut length of the path to find barcode [defacutl:5000] 
        -start_nod      :[undef | all |others] the start link nod [default:auto-find uniq edge] (manual setup: 1/2,3/4,5/6); 
        -way            :[only |all ] only produce scaffold or continue to filter the complete result [only or all]
 
USAGE


my ($barcode_stat,$nod,$edge,$region) = @ARGV;
my ($name,$outdir,$chr_type,$confict,$perl_threads,$allow_erro,$sf_set,$start_nod,$way,$start_length,$repair,$bar_edge_num,$f_path);

GetOptions(
    "n:s" => \$name,
    "i:s" => \$allow_erro,
    "len:s" => \$start_length,
    "cf:s" => \$confict,
    "fix:s" => \$repair,
    "sf:s" => \$sf_set,
    "b:s" => \$bar_edge_num,
    "fp:s" => \$f_path,
    "t:s" => \$perl_threads,
    "chr_type:s" => \$chr_type,
    "start_nod:s" => \$start_nod,
    "way:s" =>\$way,
    "o:s" => \$outdir
);

$name ||='test';
$way ||='only';
$confict ||='0';
$start_length ||='5000';
$bar_edge_num ||='1';
$f_path ||='auto';
$perl_threads ||='4';
$sf_set ||='undef';
$outdir ||=getcwd;
$start_nod ||='undef';
$repair ||='undef';
$allow_erro ||='0';

die $usage if (!$barcode_stat||!$nod||!$edge||!$chr_type);
print 'Start time:'.localtime()."\n";
if ($way eq 'only'){
    print '#Start to produce the scaffold!'."\n";
}elsif($way eq 'all'){
    print '#Start to produce the scaffold and complete path!'."\n";
}else{
    print '#Erro! Please input the corret parameters!!'."\n";
    exit;
}

##### edge set file information ###################
my %edge_cn_hash;
my %edge_index_hash;
my @uniq_edge_set;
my $j=1;
open INA,$edge or die $!;
while(<INA>){
    chomp;
    next if (/^\#/);
    my ($edg,$cn) = (split,/\s+/,$_)[0,1];
    next if ($cn == 0);
    $edge_cn_hash{$edg} = $cn;
    push @uniq_edge_set,$edg if ($cn == 1);
    my ($l,$r) = split/_/,$edg;
    $edge_index_hash{$l} = ($edg);
    my $b_edge = $r.'_'.$l;
    $edge_index_hash{$r} = ($b_edge);
    $j++;
}
close INA;
my @origin_edge; #genome link start nod 
<<<<<<< HEAD
if ($chr_type eq 'linear'){
=======
if ($chr_type eq 'liner'){
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
    @origin_edge = ($uniq_edge_set[0]);
}else{
    @origin_edge = ($uniq_edge_set[1]);
}
my $this_genome_length = sum(values %edge_cn_hash); #genome length 
##############^^Step 1 stat the basic file information #############
##### nod set file information ###################
open INB,$nod or die $!;
my %nod_all;
my %nod_hash_set;
my @total_nod;
while(<INB>){
    chomp;
    my ($a_n,$num) = split/\s+/,$_;
    push @total_nod,$a_n;
    $nod_hash_set{$a_n} = '1';
    my ($le,$ri) = split/\//,$a_n;
    push @{$nod_all{$le}},($a_n);
    my $ri_an = $ri.'/'.$le;
    push @{$nod_all{$ri}},($ri_an);
}
close INB;

##### region file infomation ######################
my %region_hash;
open INC,$region or die $!;
while(<INC>){
    chomp;
    next if (/^\#/);
    next unless($_);
    my ($a_edge,$a_len) = (split/\s+/,$_)[0,2];
    $region_hash{$a_edge} = $a_len; 
}
close INC;

##### produce the all link-result ###################
print '#Next start to get the all genome link path set! loading........'."\n";
<<<<<<< HEAD
my @all_solution;
unless (-s "$outdir/$name.info"){
    open OUTC,">$outdir/$name.info" or die $!;
=======
my %dot_graph;
my @all_solution;
open OUTC,">$outdir/$name.info" or die $!;
unless (-s "$outdir/$name.info"){
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
    while(1){
        if (@origin_edge >= 5000){
            last;
        }
        my $edge1 = shift @origin_edge;
        my @this_solution_length = split/\//,$edge1;
<<<<<<< HEAD
        if (scalar(@this_solution_length) > int($this_genome_length) - $confict){
=======
        if (scalar(@this_solution_length) > int($this_genome_length/2)){
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
            my $edg_c = $this_genome_length - scalar(@this_solution_length);
            my $nod_c = control_confict($edge1,\%nod_hash_set);
            if ($nod_c + $edg_c <= $confict){
                push @all_solution,$edge1;
                print OUTC "$edge1\t$edg_c\t$nod_c\n";
<<<<<<< HEAD
=======
                if ($chr_type eq 'cycle'){
                    my $line = ("  \"$this_solution_length[-1]\"->\"$this_solution_length[0]\"");
                    $dot_graph{$line} = '1';
                }
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
            }
        }
        my $edge_rindex = (split/_/,$edge1)[-1];
        if (exists $nod_all{$edge_rindex}){
            my @this_nod = @{$nod_all{$edge_rindex}};
            foreach my $vert (@this_nod){
                my $new_nod = (split/\//,$vert)[1];
                if (exists $edge_index_hash{$new_nod}){
                    my $new_path = $edge1.'/'.$edge_index_hash{$new_nod};
                    if (control_the_cnv($new_path,\%edge_cn_hash) eq 'T'){
<<<<<<< HEAD
=======
                        my $line = ("  \"$this_solution_length[-1]\"->\"$edge_index_hash{$new_nod}\"\[label=\"$vert\"\]");
                        $dot_graph{$line} = '1';   
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
                        push @origin_edge,$new_path;
                    }
                }
            }
        }
        if (@origin_edge == 0){
            last;
        }
    }
<<<<<<< HEAD
    ### produce more than 5000 solution result 
    if (@origin_edge < 5000){
        print 'The solution less than 5000, end !!!'."\n";
=======
    #### produce the dot graph 
    open OUTD,">$outdir/$name.2.dot" or die $!;
    print OUTD "digraph $name\_dot_graph{\n";
    print OUTD '  rankdir=LR'."\n";
    if ($chr_type eq 'cycle'){
        print OUTD '  '."\"$uniq_edge_set[1]\"".'[color=red]'."\n";
    }else{
        print OUTD '  '."\"$uniq_edge_set[0]\"".'[color=red]'."\n";
        print OUTD '  '."\"$uniq_edge_set[-1]\"".'[color=red]'."\n";
    }
    ### produce more than 5000 solution result 
    if (@origin_edge < 5000){
        print 'The solution less than 5000, end !!!'."\n";
        print OUTD join("\n",keys %dot_graph)."\n";
        # print OUTD 'labeloc=lt'."\n";
        print OUTD "label = \"Chromosome length: $this_genome_length\"\n";
        print OUTD '}'."\n";
        close OUTD;
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
    }else{
        # four threads prepare
        print 'The solution more than 5000 and it will use '.$perl_threads.' threads to run !'."\n";
        my $element_number = scalar(@origin_edge);
        for (my $i=0;$i<$element_number;$i=$i+$perl_threads){
            my $end = $i+$perl_threads;
            if ($end > $element_number){
                $end = $element_number;
            }
            my @this_some_path = @origin_edge[$i..$end-1];
            my $s_num = scalar(@this_some_path);
            my @threads_set;
            for(my $t=0;$t<$s_num;$t++){
                my ($thr) = threads->create(\&right_direct_link,$this_some_path[$t],\%edge_index_hash,\%nod_all,\%edge_cn_hash,\%nod_hash_set,$confict);
                push @threads_set,$thr;
            }
            foreach my $result (@threads_set){
                my @returnData = $result->join();
                my $this_t_num = scalar(@returnData);
                #print 'This threads produce solution number: '.$this_t_num."\n";
                push @all_solution,@returnData;
            }
        }
    }
    if (@all_solution == 0){
        print 'ERROR: Not have solution produce !!!'."\n";
        print 'Please to check the edge and nod file !!! exit'."\n";
        exit;
    }
}else{
    open INE,"$outdir/$name.info" or die $!;
    while(<INE>){
        chomp;
        if ($_){
            push @all_solution,$_;
        }
    }
    close INE;
    if (@all_solution == 0){
        print 'ERROR: Not have solution produce !!!'."\n";
        print 'Please to check the info file !!! exit'."\n";
        exit;
    }
}

close OUTC;
##### stat the useful barcode number ################
my $useful_barcode_num = int(`cat $barcode_stat |wc -l`);
if ($f_path eq 'auto'){
    if ($useful_barcode_num <= 100){
        $f_path = 1;
    }elsif ($useful_barcode_num <= 500){
        $f_path = 1;
        for (my $i=100;$i<$useful_barcode_num;$i=$i+100){
            $f_path++;
        }
    }elsif($useful_barcode_num <= 1000){
        $f_path = 5;
        for (my $i=500;$i<$useful_barcode_num;$i=$i+150){
            $f_path++;
        }
    }elsif($useful_barcode_num <=2000){
        $f_path = 8;
        for (my $i=1000;$i<$useful_barcode_num;$i=$i+200){
            $f_path++;
        }
    }elsif($useful_barcode_num <=3000){
        $f_path = 12;
        for(my $i=2000;$i<$useful_barcode_num;$i=$i+300){
            $f_path++;
        }
    }else{
        $f_path = 20;
    }
}
print '#We use barcode support number: '.$f_path."\n";

##### the duplication type check ####################
my %cn_dupli_type = check_duplication_type(\%edge_cn_hash,\%nod_all,\%nod_hash_set);
if ($repair){
    my @type_set = split/\;/,$repair;
    foreach my $ele (@type_set){
        my @the_set = split/\,/,$ele;
        my $head = shift @the_set;
        @{$cn_dupli_type{$head}} = (@the_set);
    }
}
print 'This genome contig duplication case: '."\n";
print_hash(\%cn_dupli_type,'ar');
################^^Step 2 get the start set ############
############### get the start nod ###################
my @c_set;
foreach my $t_edge (@uniq_edge_set){
    my ($ldex,$rdex) = split/_/,$t_edge;
    my $left_edge;
    my $right_edge;
    if (exists $nod_all{$ldex}){
        my @left_set = @{$nod_all{$ldex}};
        my $true_left = (split/\//,$left_set[0])[1];
        $left_edge = judge_odd_even($true_left,'l');
    }
    if (exists $nod_all{$rdex}){
        my @right_set = @{$nod_all{$rdex}};
        my $true_right = (split/\//,$right_set[0])[1];
        $right_edge = judge_odd_even($true_right,'r');
    }
    
    if ($left_edge and !$right_edge){
        my $this_path = $left_edge.'/'.$t_edge;
        my $rever_this_path = link_reverse($this_path);
        push @c_set,$rever_this_path;
    }elsif(!$left_edge and $right_edge){
        my $this_path = $t_edge.'/'.$right_edge;
        push @c_set,$this_path;
    }elsif($left_edge and $right_edge){
        my $this_path = $left_edge.'/'.$t_edge.'/'.$right_edge;
        push @c_set,$this_path;
        # my $rever_this_path = link_reverse($this_path);
        # push @c_set,($this_path,$rever_this_path);
    }else{
        print '!left and right_edge is not exists! This is not uniq edge!'."\n";
    }
}
my @rmdup = rmduplicate(\@c_set);
my @this_f_set;
my $geshu=2;
while(1){
    my %line_result;
    foreach my $this_one (@rmdup){
        my $l_num=0;
        foreach my $this_two (@rmdup){
            next if ($this_one eq $this_two);
            my $l_s = link_headtail($this_one,$this_two,$geshu);
            if ($l_s ne 'undef'){
                $line_result{$l_s} = '1';
                $l_num++;
            }
        }
        if ($l_num == 0){
            push @this_f_set,$this_one;
        }
    }
    if (keys %line_result == 0){
        last;
    }else{
        @rmdup =keys %line_result;
    }
    $geshu++;
}

my @rm_this_f_set = rmduplicate(\@this_f_set);

## reverse the direction
my @rep_set; 
foreach my $s_set (@rm_this_f_set){
    my $h_point = (split/_/,$s_set)[0];
    if (exists $nod_all{$h_point}){
        my $fan_solution = link_reverse($s_set);
        push @rep_set,($s_set,$fan_solution);
    }else{
        push @rep_set,$s_set;
    }
}
# p(@rep_set);exit;
############ complete the nod path ###############
my %path_set;
if ($start_nod eq 'undef'){
    if (@rep_set){
        foreach my $line (@rep_set){
            $path_set{$line} = '1';
        }
    }else{
        print '!Not exists uniq edge Please order the nod! You can manual setup the start_nod set! exit!!!'."\n";
        exit;
    }
}else{
    my @alon_nod;
    if ($start_nod eq 'all'){
        @alon_nod = @total_nod;
    }else{
        @alon_nod = split/\,/,$start_nod;
    } 
    foreach my $anod (@alon_nod){
        my @index = split/\//,$anod;
        my $left_link = judge_odd_even($index[0],'l');    
        my $right_link = judge_odd_even($index[1],'r');
        my $path = $left_link.'/'.$right_link;
        my $rever = link_reverse($path);
        $path_set{$path} = '1';
        $path_set{$rever} = '1';
    }
}

#^^Step 3  produce the useful information ################
############# use barcode information to produce scaffold ###############
CL:my @final_set = ();
my %sub_path = %path_set;
my $k_cycle=0;
if ($sf_set ne 'undef'){
    open INT,$sf_set or die $!;
    while(<INT>){
        chomp;
        push @final_set,$_;
    }
    close INT;
    print 'Next, We will use the fix scaffold as input !!! '."\n";
    goto LT;
}
print '#we use contig length: '.$start_length.' scaffold erro: '.$allow_erro."\n";
########### make link by the special barcode ############
while(1){
    my @list =();
    @list = keys %sub_path;
    foreach my $first (@list){
        ###1.This line basic information
        my @this_edge_set = re_sort_edge($first);
        my @this_nod_set = re_sort_nod($first);
        my $this_cn_number = $edge_cn_hash{$this_edge_set[-1]};
        my $this_uniq_edge_number = uniq(@this_edge_set);
        my @end_two_edge;
        my @last_nod;
        my ($head_point,$last_point) = (split/_/,$first)[0,-1];
        ### set the start edge number by the edge length ####
        if ($this_uniq_edge_number == 2){
            @end_two_edge = @this_edge_set[-2,-1];
            @last_nod = ($this_nod_set[-1]);
        }else{
            my $check_uniq_type = checkif_uniq(\@this_edge_set,\%edge_cn_hash);
            if ($check_uniq_type eq 'T'){
                for (my $p=2;$p<$this_uniq_edge_number;$p++){
                    my @to_stat_length = @this_edge_set[-$p..-1];
                    my $sub_check = checkif_uniq(\@to_stat_length,\%edge_cn_hash);
                    if ($sub_check eq 'T'){
                        my $two_edge_length = stat_length(\@to_stat_length,\%region_hash);
                        if ($two_edge_length < $start_length){
                            my @cut_set = cut_order_length(\@this_edge_set,$p+1);
                            if (@cut_set){
                                my $nodnum = pop @cut_set;
                                @end_two_edge = @cut_set;
                                @last_nod = tail $nodnum-1,@this_nod_set;
                            }else{
                                @end_two_edge = @this_edge_set[-$p..-1];
                                @last_nod = ($this_nod_set[-1]);
                            }
                        }else{
                            @end_two_edge = @this_edge_set[-$p..-1];
                            @last_nod = ($this_nod_set[-1]);
                        }
                    }else{
                        next;
                    }
                }
            }else{
                my @cut_set = cut_order_length(\@this_edge_set,3);
                if (@cut_set){
                    my $nodnum = pop @cut_set;
                    @end_two_edge = @cut_set;
                    @last_nod = tail $nodnum-1,@this_nod_set;
                }else{
                    @end_two_edge = @this_edge_set[-2,-1];
                    @last_nod = ($this_nod_set[-1]);
                }
            }
        }
        
        ###check if have enough length ########### 
        if (@this_edge_set == $this_genome_length){
            push @final_set,$first;
            delete $sub_path{$first};
            print '#This path genome length is enough! delete and output this path:'.$first."\n";
            next;
        }
        ###1.0 check if extists uniq link way ###############
        my $this_new_path = check_path_nod($first,$last_point,\%nod_all);
        #p($this_new_path);exit;
        if ($this_new_path ne 'undef'){
            if ($this_new_path eq 'BK'){
                push @final_set,$first;
                delete $sub_path{$first};
                print '#1.0 This point is the head point and tail point. End the link!'."\n";
                next;
            }
            if (control_the_cnv($this_new_path,\%edge_cn_hash) eq 'T'){
                $sub_path{$this_new_path} = '1'; 
                delete $sub_path{$first};
                print '#1.0 This path have a only link-nod!'."\n";
                print '#This new path is '.$this_new_path."\n";
                next;
            }else{
                push @final_set,$first;
                delete $sub_path{$first};
                print '!This path not fit cnv information! delete this path: '.$this_new_path."\n\n";
                next;
            }
            
        }
        ###2.confirm the edge duplication type
        print '#This sub_path go to confirm the edge duplication type!'.$first."\n"; 
        if (exists $cn_dupli_type{$this_edge_set[-1]}){
            my @dupli_type_set = @{$cn_dupli_type{$this_edge_set[-1]}};
            # p(@dupli_type_set);exit;
            ####2.1
            if ($dupli_type_set[0] eq 'T' and $dupli_type_set[1] eq 'T'){
                ###get the 3X3 
                my %this_barcode_hash = same_class_barcode($barcode_stat,\@end_two_edge,\@this_nod_set);
                my @point_index = split/_/,$this_edge_set[-1];
                my $copy_point;
                if ($point_index[0] > $point_index[1]){
                    $copy_point = $point_index[1].'/'.$point_index[0];
                }else{
                    $copy_point = $point_index[0].'/'.$point_index[1];
                }
                my $tandem_num=0;
                my $untandem_num=0;
                my @number;
                my $bs_num = values %this_barcode_hash;
                print '#This sub_path next link is the T-T type and go on confirm this location duplication type!'."\n";
                print '#possible tandem point is '.$copy_point."\n";
                print '#get the useful 3X3 barcode number is '.$bs_num."\n";
                
                foreach my $keys (keys %this_barcode_hash){
                    my @this_barcode_set = @{$this_barcode_hash{$keys}};
                    next unless ($this_barcode_set[1] =~/$last_point/);
                    my @anod_set = split/\,/,$this_barcode_set[1];
                    my @rm_copy_set = grep {$_ ne $copy_point} @anod_set;
                    if (grep {$_ eq $copy_point} @anod_set){
                        $tandem_num++;
                        my @bili_set = split/\;/,$this_barcode_set[2];
                        my @order_edge = grep {$_ =~ /$this_edge_set[-1]/} @bili_set;
                        my $this_number = (split/\:/,$order_edge[0])[1];
                        push @number,$this_number;
                    }
                    if (grep {$_ =~ /$last_point/} @rm_copy_set){
                        $untandem_num++
                    }
                }
                print '#get the support tandem barcode number is '.$tandem_num."\n";
                print '#get the support untandem barcode number is '.$untandem_num."\n";
                ###2.1.1
                if ($tandem_num == 0 and $untandem_num == 0){
                    print '!due to lack vital barcode! This link is break!'."\n\n";
                    #print '#find barcode condition:'.'-edge '.join(",",@end_two_edge)."\t".'-nod '.$last_nod[0]."\n";
                    push @final_set,$first;
                    delete $sub_path{$first};
                    next;
                ###2.1.2
                }elsif($tandem_num == 0 and $untandem_num > 0){
                    print '#We sub_path next link is F-T'."\n";
                    my %new_barcode_hash = class_barcode($barcode_stat,\@end_two_edge,\@last_nod,$bar_edge_num);
                    my %this_bs_stat = untandem_link($first,\%new_barcode_hash,\@last_nod,$last_point,\%edge_cn_hash,$copy_point);
                    print '#We get the assembly path and support follow: '."\n";
                    print_hash(\%this_bs_stat,'scalar');
                    my @the_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                    print '#This sub_path location is untandeme!'."\n";
                    if (exists $the_link_set[0]){
                        my @uniq_set = uniq (@the_link_set);
                        if (@uniq_set == 1){
                            $sub_path{$uniq_set[0]} = '1';
                            print '#2.1.2 We get the possible link path: '.$uniq_set[0]."\n\n";
                        }else{
                            push @final_set,$first;
                            print '#Mutiple link result. link failed!!!'."\n";
                            print join(",",@uniq_set)."\n\n";
                            # print '#This barcode set:'."\n";
                            # print_hash(\%new_barcode_hash,'ar');
                            # print_hash(\%this_bs_stat,'scalar');
                        }
                    }else{
                        ###3X3
                        my %this_bs_stat = untandem_link($first,\%this_barcode_hash,\@last_nod,$last_point,\%edge_cn_hash,$copy_point);
                        my @new_the_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                        print '!Due to lack the 3X2 barcode support.We continue use 3X3 barcode support!'."\n";
                        if (exists $new_the_link_set[0]){
                            my @uniq_set = uniq (@new_the_link_set);
                            if (@uniq_set == 1){
                                $sub_path{$uniq_set[0]} = '1';
                                print '#We get the assembly path and support follow: '."\n";
                                print_hash(\%this_bs_stat,'ar');
                                print '#2.1.2 We get this cycle final link path: '.$uniq_set[0]."\n\n";
                            }else{
                                push @final_set,$first;
                                print '#Mutiple link result. link failed!!!'."\n\n";
                                print join(",",@uniq_set)."\n";
                                # print '#This barcode set:'."\n";
                                # print_hash(\%this_barcode_hash,'ar');
                                # print_hash(\%this_bs_stat,'scalar');
                            }
                        }else{
                            push @final_set,$first;
                            print '#No result! link failed!!!'."\n\n";
                            # print '#This barcode set:'."\n";
                            # print_hash(\%this_barcode_hash,'ar');
                        }
                    }
                ### 2.1.3 
                }elsif($tandem_num > 0 and $untandem_num == 0){
                    print '#This sub_path next link is T-F'."\n";
                    my $average_cn = floor(sum(@number)/scalar(@number));
                    if ($average_cn >= $this_cn_number - 1 or $average_cn < 1){
                        print '!Why this rpkm_cn canclute erro?'."\n";
                        print $this_edge_set[-1]."\t".'conclute cn:'.$average_cn."\t".'first cn:'.$this_cn_number."\n";
                        push @final_set,$first;
                        delete $sub_path{$first};
                        print '#link failed!!!'."\n";
                        next;
                    }
                    my $new_path = $first;
                    my $last_edge = (split/\//,$first)[-1];
                    for (my $j=0;$j<$average_cn-1;$j++){
                        $new_path = $new_path.'/'.$last_edge;
                    }
                    print '#This location is the tandem!'."\n";
                    print '#We get the tandem path result:'.$new_path."\n";
                    print '#next we go on the next link:'."\n";
                    my %this_bs_stat = untandem_link($new_path,\%this_barcode_hash,\@last_nod,$last_point,\%edge_cn_hash);
                    my @inthis_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                    if (exists $inthis_link_set[0]){
                        my @uniq_set = uniq (@inthis_link_set);
                        if (@uniq_set == 1){
                            $sub_path{$uniq_set[0]} = '1';
                            print '#We get the assembly path and support follow: '."\n";
                            print_hash(\%this_bs_stat,'scalar');
                            print '#2.1.3 We get this cycle final path:'.$uniq_set[0]."\n\n";
                        }else{
                            print '#Mutiple link result. link failed!!!'."\n";
                            print join(",",@uniq_set)."\n\n";
                            # print '#This barcode set:'."\n";
                            # print_hash(\%this_barcode_hash,'ar');
                            # print_hash(\%this_bs_stat,'scalar');
                            push @final_set,$new_path;
                        }
                    }else{
                        print '#No result. link failed!!!'."\n\n";
                        # print '#This barcode set:'."\n";
                        # print_hash(\%this_barcode_hash,'ar');
                        push @final_set,$new_path;
                    }  
                }else{
                    print '!exists the tandem and untandem barcode? This barcode can\'t know this lcation infroamtion!'."\n\n";
                    print '#This path is '.$first."\n";
                    print '#This barcode set:'."\n";
                    print_hash(\%this_barcode_hash,'ar');
                    delete $sub_path{$first};
                    push @final_set,$first;
                }
            ### 2.2
            }elsif($dupli_type_set[0] eq 'T' and $dupli_type_set[1] eq 'F'){
                print '#This sub_path next link is the T-F'."\n";
                my @dup_link;
                for (my $i=0;$i<$this_cn_number-1;$i++){
                    push @dup_link,(split/\//,$first)[-1];
                }
                my $dup_edge = join("/",@dup_link);
                my $link_result = $first.'/'.$dup_edge;
                print '#2.2 We get the complete tandem duplication path:'.$link_result."\n";
                my %inthis_t_bs = same_class_barcode($barcode_stat,\@end_two_edge,\@last_nod);
                unless(%inthis_t_bs){
                    print '!lack the some barcode:'.$link_result."\n";
                    print '#find barcode condition:'.'-edge '.join(",",@end_two_edge)."\t".'-nod '.join(",",@last_nod)."\n";
                    push @final_set,$link_result;
                    delete $sub_path{$first};
                }
                my %this_bs_stat = untandem_link($link_result,\%inthis_t_bs,\@last_nod,$last_point,\%edge_cn_hash);
                my @inthis_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                print '#next we go on the next link:'."\n";
                if (exists $inthis_link_set[0]){
                    my @uniq_set = uniq (@inthis_link_set);
                    if (@uniq_set == 1){
                        $sub_path{$uniq_set[0]} = '1';
                        print '#We get the assembly path and support follow: '."\n";
                        print_hash(\%this_bs_stat,'scalar');
                        print '#2.2 We get this cycle final path:'.$uniq_set[0]."\n\n";
                    }else{
                        push @final_set,$link_result;
                        print '#Mutiple link result. link failed!!!'."\n";
                        print join(",",@uniq_set)."\n\n";
                        # print '#This barcode set:'."\n";
                        # print_hash(\%inthis_t_bs,'ar');
                        # print_hash(\%this_bs_stat,'scalar');
                    }
                }else{
                    push @final_set,$link_result;
                    print '#No link result. link failed!!!'."\n\n";
                    # print '#This barcode set:'."\n";
                    # print_hash(\%inthis_t_bs,'ar');
                }
            ### 2.3
            }elsif($dupli_type_set[0] eq 'F' and $dupli_type_set[1] eq 'T'){
                print '#This sub_path next link is the F-T'."\n";
                my %this_barcode_hash = class_barcode($barcode_stat,\@end_two_edge,\@last_nod,$bar_edge_num);
                my %this_bs_stat = untandem_link($first,\%this_barcode_hash,\@last_nod,$last_point,\%edge_cn_hash);
                if ($head_point == $last_point){    
                    my $rever_first = link_reverse($first);
                    my @the_head_edge_set = re_sort_edge($rever_first);
                    my @the_head_nod_set = re_sort_nod($rever_first);
                    my $tail_num = scalar(@end_two_edge);
                    my @head_cut_set = cut_order_length(\@the_head_edge_set,$tail_num);
                    my $nnodnum = pop @head_cut_set;
                    my @head_nod = tail $nnodnum-1,@the_head_nod_set;
                    my %head_barcode_hash = class_barcode($barcode_stat,\@head_cut_set,\@head_nod,$bar_edge_num);
                    my %head_bs_stat = untandem_link($rever_first,\%head_barcode_hash,\@head_nod,$last_point,\%edge_cn_hash);
                    my %this_rl = fix_the_doublediret(\%head_bs_stat,\%this_bs_stat);
                    my @rl_return_set = keys %this_rl;
                    if (@rl_return_set == 1){
                        $sub_path{$rl_return_set[0]} = '1';
                        print '#We get the assembly path and support follow: '."\n";
                        print_hash(\%this_bs_stat,'scalar');
                        print '#2.3 We get this cycle final path:'.$rl_return_set[0]."\n\n";
                    }else{
                        push @final_set,$first;
                        print '#Mutiple link result. link failed!!!'."\n";
                        print join(",",@rl_return_set)."\n\n";
                        # print '#This barcode set:'."\n";
                        # print_hash(\%this_barcode_hash,'ar');
                        # print_hash(\%this_bs_stat,'scalar');
                    }
                }else{
                    my @the_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                    if (exists $the_link_set[0]){
                        my @uniq_set = uniq (@the_link_set);
                        if (@uniq_set == 1){
                            $sub_path{$uniq_set[0]} = '1';
                            print '#We get the assembly path and support follow: '."\n";
                            print_hash(\%this_bs_stat,'scalar');
                            print '#2.3 We get this cycle final path:'.$uniq_set[0]."\n\n";
                        }else{
                            push @final_set,$first;
                            print '#Mutiple link result. link failed!!!'."\n";
                            print join(",",@uniq_set)."\n\n";
                            # print '#This barcode set:'."\n";
                            # print_hash(\%this_barcode_hash,'ar');
                            # print_hash(\%this_bs_stat,'scalar');
                        }
                    }else{
                        push @final_set,$first;
                        print '#No link result. link failed!!!'."\n\n";
                        # print '#This barcode set:'."\n";
                        # print_hash(\%this_barcode_hash,'ar');
                    }
                }
            }else{
                print '#It\'s impossible for the duplication type!'."\n";
            }
        ### 2. the alon edge link ##########################
        }else{
            print '#This link only one-nod:'."\n";
            if (exists $nod_all{$last_point}){
                my @cnod = @{$nod_all{$last_point}};
                if (@cnod > 1){
                    print 'It is wrong!!!'."\n";
                    print join(",",@cnod)."\n";
                }
                my $rdex = (split/\//,$cnod[0])[1];
                my $next_edge = $edge_index_hash{$rdex};
                my $new_path = $first.'/'.$next_edge;
                $sub_path{$new_path} = '1';
            print '#We get new path:'.$new_path."\n\n";
            }else{
                push @final_set,$first;
                delete $sub_path{$first};
                print '#This point is the head point and tail point. End the link!'."\n";
            }
            
        }
        delete $sub_path{$first} if (exists $sub_path{$first});
    }
    $k_cycle++;
    last if (keys %sub_path == 0);
}

###### rm the duplication path and merge the common path #########
print "!!!Finally $k_cycle cycle in this!\n\n";

my @output_set;
foreach my $this_solution (@final_set){
    my $this_edge_number = split/\//,$this_solution;
    if ($this_edge_number > int($this_genome_length/2)){
        my $edge_confict = $this_genome_length - $this_edge_number;
        my $nod_confict = control_confict($this_solution,\%nod_hash_set);
        if ($edge_confict + $nod_confict <= $confict){
            push @output_set,$this_solution;
        }
    }
}
LT:if ($sf_set eq 'undef'){
    open OUTB,">$outdir/$name.result" or die $!;
    if (@output_set){
        my @uniq_output_set = rmduplicate(\@output_set);
        my $f_s_n = scalar(@uniq_output_set);
        my $a_s_n = scalar(@all_solution);
        print OUTB join("\n",@uniq_output_set)."\n";
        print '#exists complete path! This link is end.'."\n";
        print 'ALL solution number: '.$a_s_n."\t".'Filter solution number: '.$f_s_n."\n";
    }else{
        print '#Before the rmdup scaffold set: '.join("\t",@final_set)."\n";
        my @rm_final_set = rmduplicate(\@final_set);
        if (@rm_final_set){
            print '#Output the scaffold set:'."\n";
            print join("\t",@rm_final_set)."\n";
            open OUTA,">$outdir/$name.scaffold" or die $!;
            print OUTA join("\n",@rm_final_set)."\n";
            close OUTA;
        }else{
            print '!!!Not exists scaffold. You should adjust the start_nod!!!'."\n";
        }
        if ($way eq 'only'){
            exit;
        }
        my @filter_set;
        my $solu_number=0;
        foreach my $apath (@all_solution){
            my $check_result = check_solution($apath,\@rm_final_set,$allow_erro,$chr_type);
            if ($check_result eq 'T'){
                $solu_number++;
            }
        }
        if ($solu_number != 0){
            my $a_s_n = scalar(@all_solution);
            print 'ALL solution number: '.$a_s_n."\t".'Filter solution number: '.$solu_number."\n";
        }else{
            ## upadata the paramaters
            if ($allow_erro == 1){
                if ($f_path == 1){
                    print 'No result !!! It is a wrong !!! please adjust the paramaters!!!'."\n";
                    goto END;
                }else{
                    $f_path = $f_path -2;
                    if ($f_path == 0){
                        $f_path = 1;
                    }
                }
            }
            if ($start_length > 21000){
                $start_length = 5000;
                $allow_erro = 1;
            }else{
                $start_length = $start_length + 5000;
            }
            goto CL;
        }
    }
}else{
    my @fix_set;
    open OUTB,">$outdir/$name.result" or die $!;
    my $f_s_n=0;
    foreach my $apath (@all_solution){
        my $check_result = check_solution($apath,\@final_set,$allow_erro,$chr_type);
        if ($check_result eq 'T'){
            $f_s_n++;
        }
    }
    my $a_s_n = scalar(@all_solution);
    print OUTB join("\n",@fix_set)."\n";
    print 'ALL solution number: '.$a_s_n."\t".'Filter solution number: '.$f_s_n."\n";
}
END:print '****************************This work is finished!!!*************************************'."\n";
close OUTB;

print 'End time:'.localtime()."\n";
###############sub check_path_nod #############
sub check_path_nod{
    my ($sub,$last,$all_nod) = @_;
    my @this_path_nod_set = re_sort_nod($sub);
    my @end_nod_set;

    if (exists ${$all_nod}{$last}){
        @end_nod_set = @{${$all_nod}{$last}};
    }else{
        return 'BK';
        last;
    }
    my $new_path = 'undef';
    if (@end_nod_set == 1){
        my ($ldex,$rdex) = split/\//,$end_nod_set[0]; 
        my $right_edge;
        if ($last == $ldex){
            $right_edge = judge_odd_even($rdex,'r');
        }else{
            $right_edge = judge_odd_even($ldex,'r');
        }
        $new_path = $sub.'/'.$right_edge;
    }
    else{
        my @sort_nod_set = re_sort_nod(\@end_nod_set);
        my @jiao_ji = intersect(@this_path_nod_set,@sort_nod_set);
        if (@jiao_ji){
            my @left_set = array_diff(@sort_nod_set,@jiao_ji);
            if (@left_set == 1){
                my ($ldex,$rdex) = split/\//,$left_set[0]; 
                my $right_edge;
                if ($last == $ldex){
                    $right_edge = judge_odd_even($rdex,'r');
                }else{
                    $right_edge = judge_odd_even($ldex,'r');
                }
                $new_path = $sub.'/'.$right_edge;
            }   
        }
    } 
    return $new_path;
}

###############sub untandem_link ##############
sub untandem_link{
    my ($aedge,$the_barcode_data_set,$the_bkp,$endnod,$cn_hash,$cp) = @_;
    ###barcode link set #########
    my %bs_stat;
    my %support_barcode_position;
    foreach my $keys (keys %$the_barcode_data_set){
        my @the_set = @{${$the_barcode_data_set}{$keys}};
        my @the_barcode_edge_set = split/\,/,$the_set[0];
        my @the_barcode_nod_set = split/\,/,$the_set[1];
        my @diff_nod = array_diff(@the_barcode_nod_set,@$the_bkp);
        my @one_nod;
        if ($cp){
            @one_nod = grep {$_ ne $cp} @diff_nod;
        }
        @one_nod = grep {$_ =~ /$endnod/} @diff_nod;
        if (@one_nod){
            foreach my $a_nod (@one_nod){
                my $link_result = once_link($aedge,\@the_barcode_edge_set,$a_nod);
                if ($link_result ne 'undef'){
                    if (control_the_cnv($link_result,\%$cn_hash) eq 'T'){
                        if (exists $bs_stat{$link_result}){
                            $bs_stat{$link_result}++;
                        }else{
                            $bs_stat{$link_result} = 1;
                        }
                        push @{$support_barcode_position{$link_result}},$keys;
                    }
                }
            }
        }
    }
    foreach my $doublekey (keys %support_barcode_position){
        my @bid = @{$support_barcode_position{$doublekey}};
        print "$doublekey\t".join(",",@bid)."\n";
    }
    return %bs_stat;
}
#############sub right_direct_link ################ 
sub right_direct_link{
    my ($origin,$edge_set,$nod_set,$edge_hash,$anodset,$allow_confict)= @_;
    my $genome_length = sum(values %$edge_hash);
    my @solution;
    my @list = ($origin);
    while(1){
        my $edge1 = shift @list;
        my $this_solution_length = split/\//,$edge1;
        if ($this_solution_length > int($genome_length/2)){
            my $edg_c = $genome_length-$this_solution_length;
            my $nod_c = control_confict($edge1,\%$anodset);
            if ($edg_c+$nod_c <= $allow_confict){
                push @solution,$edge1;
                print OUTC "$edge1\t$edg_c\t$nod_c\n";
            }
        }
        my $edge_rindex = (split/_/,$edge1)[-1];
        if (exists ${$nod_set}{$edge_rindex}){
            my @this_nod = @{${$nod_set}{$edge_rindex}};
            foreach my $vert (@this_nod){
                my $new_nod = (split/\//,$vert)[1];
                if (exists ${$edge_set}{$new_nod}){
                    my $new_path = $edge1.'/'.${$edge_set}{$new_nod};
                    if (control_the_cnv($new_path,\%$edge_hash) eq 'T'){
                            push @list,$new_path;
                    }
                }
            }
        }
        if (@list == 0){
            last;
        }
    }
    return @solution;
}

###############sub rmduplicate ###################
sub rmduplicate{
    my $input = $_[0];
    my @uniq = uniq(@$input);
    my %hash_set;
    foreach my $sub (@uniq){
        $hash_set{$sub} = '1';
    }
    foreach my $lu1 (@uniq){
        my $len1 = split/\//,$lu1;
        foreach my $lu2 (@uniq){
            next if ($lu1 eq $lu2);
            my $len2 = split/\//,$lu2;
            if ($len1 < $len2){
                if ($lu2 =~ /$lu1/){
                    delete $hash_set{$lu1};
                }else{
                    my $rever_lu1 = link_reverse($lu1);
                    if ($lu2 =~ /$rever_lu1/){
                        delete $hash_set{$lu1};
                    }
                }
            }elsif($len1 > $len2){
                if ($lu1 =~ /$lu2/){
                    delete $hash_set{$lu2};
                }else{
                    my $rever_lu2 = link_reverse($lu2);
                    if ($lu2 =~ /$rever_lu2/){
                        delete $hash_set{$lu2};
                    }
                }
            }else{
                my $rever_lu2 = link_reverse($lu2);
                if ($lu1 =~ /$rever_lu2/){
                    if (exists $hash_set{$lu1} and exists $hash_set{$lu2}){
                        delete $hash_set{$lu2};
                    }
                }
            }
        }
    }
    my @output = keys %hash_set;
    my @filter_set;
    foreach my $element (@output){
        my $edge_number = split/\//,$element;
        if ($edge_number > 2){
            push @filter_set,$element;
        }
    }
    return @filter_set;
}

######################sub 
sub once_link{
    my ($origin,$aedge_set,$anod)=@_;
    my $new_path = 'undef';      
    my $edge_rindex = (split/_/,$origin)[-1];
    my ($nod_lindex,$nod_rindex)=split/\//,$anod;
    my $new_nod;
    if ($edge_rindex == $nod_lindex or $edge_rindex == $nod_rindex){
        $new_nod = $nod_rindex if ($edge_rindex == $nod_lindex);
        $new_nod = $nod_lindex if ($edge_rindex == $nod_rindex);
        foreach my $aedge (@$aedge_set){
            my @index2=split/_/,$aedge;
            my $lindex =$index2[0];
            my $rindex = $index2[-1];
            my $new_edge;                      
            if ($lindex == $new_nod or $rindex == $new_nod){
                $new_edge = $aedge if ($lindex == $new_nod);
                $new_edge = link_reverse($aedge) if ($rindex == $new_nod);
                $new_path = $origin.'/'.$new_edge;
            }
        }
    }
    return $new_path;
}

########re-sort the sub_path edge set likes 11_10→10_11 ##########
sub re_sort_edge{
    my $input = $_[0];
    my @array = split/\//,$input;
    my @new_set;
    foreach my $edge (@array){
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

##################################################################
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

##########sub control_the_cnv($new_path,\%edge) ####################
sub control_the_cnv{
    my ($this_path,$this_hash)=@_;
    my @set = re_sort_edge($this_path);
    my @uniq =uniq(@set);
    my @stat;
    my $check = 'T';
    foreach my $re_ed (@uniq){
        my $num = grep {$_ eq $re_ed} @set;
        if ($num <= ${$this_hash}{$re_ed}){
            next; 
        }else{
            $check = 'F';
            last;
        }
    }
    return $check;
}

###############sub link_reverse ################
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

############################################
sub check_duplication_type{
    my ($edge_hash,$nod_all_hash,$nod)=@_;
    my %nod_hash = %$nod_all_hash;
    my @hight_duplica_set;
    my %output;
    foreach my $keys (keys %$edge_hash){
        my $num = ${$edge_hash}{$keys};
        next if ($num < 2);
        my @aset = split/_/,$keys;
        my $this_nod = $aset[0].'/'.$aset[1];
        my $tandem = 'F';
        my $untandem = 'F'; 
        if ($num == 2){
            if (exists ${$nod}{$this_nod}){
                $tandem = 'T';
                $untandem = 'F';
            }else{
                $untandem = 'T';
                $tandem = 'F';
            }
        }elsif($num > 2){
            my @left = @{$nod_hash{$aset[0]}};
            my $left_num = scalar(@left);
            my @right = @{$nod_hash{$aset[1]}};
            my $right_num = scalar(@right);
            my $common = 0;
            if (exists ${$nod}{$this_nod}){
                $common = 1;
            }
            if ($common == 1){
                if ($left_num == 2 and $right_num ==2){
                    $tandem = 'T';
                    $untandem = 'F';
                }elsif($left_num >2 or $right_num >2){
                    $tandem = 'T';
                    $untandem = 'T';
                }
            }else{
                $tandem = 'F';
                $untandem = 'T';
            }
        }
        @{$output{$keys}} = ($tandem,$untandem);
    }
    return %output;
}

###############################################################
sub print_hash{
    my ($input,$at) = @_;
    if (%$input){
        foreach my $keys (keys %$input){
            if ($at eq 'ar'){
                my @this_set = @{${$input}{$keys}};
                print "$keys\t".join("\t",@this_set)."\n";
            }else{
                my $this_element = ${$input}{$keys};
                print "$keys\t$this_element\n";
            }
        }
    }else{
        print 'Empty!'."\n";
    }
    print "\n";
}

#################sub control_confict #######################
sub control_confict{
    my ($edg,$nod_set) = @_;
    my @nodset = split/_/,$edg;
    my $first_poin = shift @nodset; 
    my $last_poin = pop @nodset;
    my $local;
    if ($first_poin > $last_poin){
        $local = $last_poin.'/'.$first_poin;
    }else{
        $local = $first_poin.'/'.$last_poin;
    }
    my %double_nod = %$nod_set;
    delete $double_nod{$local} if (exists $double_nod{$local});
    foreach my $nd (@nodset){
        my($ln,$rn)=split/\//,$nd;
        if ($ln > $rn){
            $nd = $rn."/".$ln;
        }
        if ($double_nod{$nd}){
            delete $double_nod{$nd} if (exists $double_nod{$nd});
        }
    }
    my $nod_confict= values %double_nod;
    return $nod_confict;
}

################################################################
sub cut_order_length{
    my ($input,$wtnumber) = @_;
    my $edge_num = scalar(@$input);
    my $uniq_path_number = scalar(uniq(@$input)); 
    my @ct_set;
    if ($wtnumber > $uniq_path_number){
        return @ct_set;
        last;
    }
    my $w;
    for ($w=$wtnumber;$w<=$edge_num;$w++){
        my @cut_set = tail $w,@$input;
        my @uniq_cut_set = uniq(@cut_set);
        if (@uniq_cut_set == $wtnumber){
            @ct_set = @uniq_cut_set;
            last;
        }
    }
    push @ct_set,$w;
    return @ct_set;
}

####################################################
sub stat_length{
    my ($apath,$length_hash) = @_;
    my $the_length = 0;
    foreach my $t_path (@$apath){
        $the_length = $the_length + ${$length_hash}{$t_path};
    }
    return $the_length;
}

############################################################
sub barcode_support_num{
    my ($input,$filter) = @_;
    my @return_set;
    foreach my $keys (keys %$input){
        my $this_set_num = ${$input}{$keys};
        if ($this_set_num < $filter){
            next;
        }else{
            push @return_set,$keys;
        }
    }
    return @return_set;
}

####################################################
sub check_solution{
    my ($path,$the_scaf,$supp_num,$chrmo_type) = @_;
    my $scaf_num = scalar(@$the_scaf);
    my $stat_num=0;
<<<<<<< HEAD
    if ($chrmo_type eq 'linear'){
=======
    if ($chrmo_type eq 'liner'){
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
        foreach my $this_local (@$the_scaf){
            if ($path =~ /$this_local/){
                $stat_num++;
            }else{
                my $rever_local = link_reverse($this_local);
                if ($path =~ /$rever_local/){
                    $stat_num++;
                }else{
                    next;
                }
            }
        }
    }elsif ($chrmo_type eq 'cycle'){
        foreach my $this_local (@$the_scaf){
            if ($path =~ /$this_local/){
                $stat_num++;
            }else{
                my $rever_edge = link_reverse($this_local);
                if ($path =~ /$rever_edge/){
                    $stat_num++;
                }else{
                    my ($tail_point,$head_point) = (split/_/,$path)[0,-1];
                    my $this_nod = $head_point.'/'.$tail_point;
                    if ($this_local =~ /$this_nod/){
                        my ($head,$tail) = split/$this_nod/,$this_local,2;
                        my $new_head = $head.$head_point;
                        my $new_tail = $tail_point.$tail;
                        if ($path =~ /^$new_tail/ and $path =~/$new_head$/){
                            $stat_num++;
                        }
                    }else{
                        if ($rever_edge =~ /$this_nod/){
                            my ($head,$tail) = split/$this_nod/,$rever_edge,2;
                            my $new_head = $head.$head_point;
                            my $new_tail = $tail_point.$tail;
                            if ($path =~ /^$new_tail/ and $path =~/$new_head$/){
                                $stat_num++;
                            }
                        }else{
                            next;
                        }
                    }
                }
            }
        }
    }else{
        print 'please input the corret chrmosome type !!!'."\n";
        exit;
    }
    if ($stat_num >= $scaf_num - $supp_num){
        print OUTB $path."\n";
        return 'T';
    }else{
        return 'F';
    }
}

<<<<<<< HEAD
=======

>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
################################################################
sub checkif_uniq{
    my ($input,$ed_hash) = @_;
    my $ct = 'F';
    foreach my $ad (@$input){
        if (${$ed_hash}{$ad} == 1){
            $ct = 'T';
        }
    }
    return $ct;
}

<<<<<<< HEAD
=======

>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
##############################################################
sub fix_the_doublediret{
    my ($head_hash,$tail_hash) = @_;
    my @head_solution = keys %$head_hash;
    my @tail_solution = keys %$tail_hash;
    my @this_last_edge;
    foreach my $ahl (@head_solution){
        my @lls = split/\//,$ahl;
        push @this_last_edge,$lls[-1];
    }
    foreach my $asl (@tail_solution){
        my @tls = split/\//,$asl;
        if (grep {$_ eq $tls[-1]} @this_last_edge){
            delete ${$tail_hash}{$asl};
        }
    }
    return %$tail_hash;
}

<<<<<<< HEAD
# my $first_set = '2_3/4_5/6_7';
# my $second_set = '4_5/6_7/8_9/10_11';
# my $myi = 2;
# my $output = link_headtail($first_set,$second_set,$myi);
## output = 2_3/4_5/6_7/8_9/10_11
=======

>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
#############sub link_headtail ###################3
sub link_headtail{
    my ($first,$second,$i) = @_;
    my @set1 = split/\//,$first;
    my $num1 = scalar(@set1);
    my @set2 = split/\//,$second;
    my $num2 = scalar(@set2);
    my $range;
    ##confirm the overlap ####
    if ($num1<$num2){
        $range=$num1;
    }else{
        $range=$num2;
    }
    my $new_path = 'undef';
    ## first ##
    my @t1 = tail $i,@set1;
    ## second ##
    my @h2 = head $i,@set2;
    ###
<<<<<<< HEAD
    my $cresult1 = compare_two_array(\@t1,\@h2);
    if ($cresult1  == 1){
=======
    if (@t1 ~~ @h2){
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
        my $p1 = ($num1-$i);
        my @link1 = head $p1,@set1;
        my @new = (@link1,@set2);
        $new_path = join("/",@new);
    }else{
        my $p1 = ($num1-$i);
        my $third = link_reverse($second);
        my @set3 = split/\//,$third;
        my @h3 = head $i,@set3;
<<<<<<< HEAD
        my $cresult2 = compare_two_array(\@t1,\@h3);
        if($cresult2 == 1){
=======
        if(@t1 ~~ @h3){
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
            my @link3 = head $p1,@set1;
            my @new = (@link3,@set3);
            $new_path = join("/",@new);
        }
        
    }
    return $new_path ;
}
<<<<<<< HEAD
####################################################
sub compare_two_array{
    my ($array1,$array2) = @_;
    my $flag = 0;
    my %hash1 = map { $_ => 1 } @$array1;
    my %hash2 = map { $_ => 1 } @$array2;
    if (keys %hash1 == keys %hash2) {
        $flag = 1;
        foreach my $key (keys %hash1) {
            unless (exists $hash2{$key}) {
                $flag = 0;
                last;
            }
        }
    }
    return $flag;
}
=======
>>>>>>> 4995c299facc3f5a9b96a35ade69e745dbe83e6a
