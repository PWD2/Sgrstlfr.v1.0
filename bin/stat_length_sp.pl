#!/usr/bin/perl -w

use strict;
use Cwd;
use threads;
use Array::Utils qw(:all);
use List::Util qw/max min/;
use List::MoreUtils qw/indexes uniq /;
use Data::Printer;
use Getopt::Long;

my ($mutiple_solution,$info,$total_set,$barcode_set) = @ARGV;

my ($minnum,$maxnum,$perl_threads,$gap_length,$chrtype,$precen_max,$sp_minus,$mbdup,$name,$outdir);
my $usage=<<USAGE;

    This script is order to decline the solution of SCRaMbLE by the appoint barcode type!!!
    v1.0 2023.2.21 pangwending
    perl $0 <result> <info> <total> <barcode(merge)> [-option]

    -n              : the prefix name [default:test]
    -o              : the outdir [default:./]
    -t              : the threads number [default:4]
    -min            : the minnum of the edge [default:3]
    -max            : the maxnum of the edge [default:12]
    -gap_len        : the gap length [default:500]
    -chr_type       : the chrmosome type [default:liner]
    -max_b_len      : the max barcode length precent of genome_length [default:0.5]
    -max_b_dup      : the max edge duplication number in one barcode [default:5]
    -sp_minus       : [max |others]the barcode support decline precent of two solution [default:0.01]

USAGE

GetOptions(
    "n:s" => \$name,
    "o:s" => \$outdir,
    "t:s" => \$perl_threads,
    "min:s" => \$minnum,
    "max:s" => \$maxnum,
    "gap_len:s" => \$gap_length,
    "chr_type:s" => \$chrtype,  
    "max_b_len:s" => \$precen_max,
    "max_b_dup:s" => \$mbdup,
    "sp_minus:s" => \$sp_minus
);

die $usage if (!$mutiple_solution||!$barcode_set||!$total_set);

$name ||='test';
$outdir ||=getcwd;
$perl_threads ||='4';
$minnum ||=3;
$maxnum ||=12;
$gap_length ||=500;
$chrtype ||='liner';
$precen_max ||=0.5;
$mbdup ||=5;
$sp_minus ||=0.01;

unless (-s $mutiple_solution){
    print '!!!The mutiple result is Empty!'."\n";
    print 'It will use the info file to filter !'."\n";
    unless (-s $info){
        print 'But the info file is not exists !!! End !'."\n";
        exit;
    }
}

open OUTB,">$outdir/$name.bs.log" or die $!;
print OUTB 'Start time: '.localtime()."\n";
print OUTB '########################### fiter the mutiple solution by the co-barcode reads will start #################################'."\n";

my @solution_set;
my $max_edge;
my $path_num=0;
if (-s $mutiple_solution){
    open INA,$mutiple_solution or die $!;
}else{
    open INA,$info or die $!;
}
while(<INA>){
    chomp;
    my $solu = (split/\s+/,$_)[0];
    push @solution_set,$solu;
    my $edge_number = split/\//,$solu; 
    $max_edge = $edge_number-1;
    $path_num++;
}
close INA;
if (@solution_set == 0){
    print 'ERROR: Not have the solution as input !!!!'."\n";
    exit;
}else{
    print OUTB 'We detact the solution number: '.$path_num."\n";
}
# p(%solution_edge_set);p(%solution_nod_set);exit;

my %region_hash;
my $genome_length;
my $first_edge;
my $last_edge;
my @short_length_edge;
open INB,$total_set or die $!;
while(<INB>){
    chomp;
    if ($_ =~ /genome/){
        $genome_length = (split/\:/,$_)[1];
    }
    next if (/^\#/);
    next unless($_);
    my ($a_edge,$a_len) = (split/\s+/,$_)[0,2];
    if ($a_len < $gap_length){
        push @short_length_edge,$a_edge;
    }
    if ($.== 2){
        $first_edge = $a_edge;
    }
    $region_hash{$a_edge} = $a_len;
    $last_edge = $a_edge; 
}
close INB;
my $max_allow_length = int($genome_length*$precen_max);
# p($first_edge);p($last_edge);p($max_allow_length);p(%region_hash);exit;


print OUTB '!we use the max barcode length: '.$max_allow_length."\t".'max gap length: (10)'.$gap_length."\n";
print OUTB '!This chrmosome type: '.$chrtype."\n";
print OUTB '!We start to use the barcode type: '.(split/\//,$barcode_set)[-1]."\n\n";

my $useful_barcode=0;
#open OUTD,">$outdir/$name.gap.distrubution.txt" or die $!;
my %bs_set;
my %cluster_set;
open INC,$barcode_set or die $!;
<INC>;
while(<INC>){
    chomp;
    my ($barcode,$edge_num,$edge_nod,$nod_num,$node,$bili) = split/\s+/,$_;
    next if ($nod_num <1);
    if ($edge_num >=$minnum and $edge_num <= $maxnum){
        $useful_barcode++;
        push @{$cluster_set{$edge_nod}{$node}},($barcode);
            
    }
}
my $uniq_barcode=0;
foreach my $keys (keys %cluster_set){
    my %second_hash = %{$cluster_set{$keys}};
    my @edge_set = split/\,/,$keys;
    foreach my $skeys (keys %second_hash){
        my @nod_set = split/\,/,$skeys;
        my @bar_set = @{$second_hash{$skeys}};
        my $b_num = scalar(@bar_set);
        $uniq_barcode++;
        @{$bs_set{$bar_set[0]}} = ([@edge_set],[@nod_set],$b_num);
    }
}

print OUTB '!The useful barcode number: '.$useful_barcode."\n";
print OUTB '!The uniq barcode type number: '.$uniq_barcode."\n";
close INB;

open OUTA,">$outdir/$name.path.log" or die $!;
my @all_set;
my @arr;
my $element_number = scalar(@solution_set);
for (my $i=0;$i<$element_number;$i=$i+$perl_threads){
    my $end = $i+$perl_threads;
    if ($end > $element_number){
        $end = $element_number;
    }
    my @this_some_path = @solution_set[$i..$end-1];
    my $s_num = scalar(@this_some_path);
    my @threads_set;
    for(my $t=0;$t<$s_num;$t++){
        my $thr = threads->create(\&stat_bs_number,$this_some_path[$t]);
        push @threads_set,$thr;
    }
    foreach my $result (@threads_set){
        $result->join();
    }
}
close OUTA;
my $second;
my $output_number=0;
open OUTD,">$outdir/$name.path" or die $!;
open IND,"sort -k 1 -n -r $outdir/$name.path.log |" or die $!;
while(<IND>){
    chomp;
    my ($bs,$solution) = split/\s+/,$_;
    if ($.==1){
        $second = $bs;
    }
    my $reduce = $second - $bs;
    my $bili = $reduce/$second;
    if ($sp_minus eq 'max'){
        if ($bili == 0){
            print OUTD "$bs\t$solution\n";
        }else{
            last;
        }
    }else{
        if ($bili <= $sp_minus and $second > 0){
            print OUTD "$bs\t$solution\n";
            $output_number++;
        }else{
            last;
        }
    }
    $second = $bs;
}
close OUTD;
close IND;
print OUTB "\n".'!We get the solution number: '.$output_number."\n";
print OUTB 'Start time: '.localtime()."\n";
print OUTB 'End time:'.localtime()."\n";



close OUTB;


###########################################################################################################
sub stat_bs_number{
    my $solution = $_[0];
    my $bs_stat=0;
    my @path_set = re_sort_edge($solution);
    my @nd_set = re_sort_nod($solution);
    foreach my $barcode (keys %bs_set){
        ### the basic solution index
        my @edge_set = @{${$bs_set{$barcode}}[0]};
        my @nod_set = @{${$bs_set{$barcode}}[1]};
        
        my %this_edge_index_hash = re_index_edge(\@path_set,\@edge_set);
        my @this_nod_index_set = re_index_nod(\@nd_set,\@nod_set);
        
        ### uniq edge index add
        foreach my $mkey (keys %this_edge_index_hash){
            my @the_s_num = @{$this_edge_index_hash{$mkey}};
            if (@the_s_num == 1){
                push @this_nod_index_set,$the_s_num[0];
                delete $this_edge_index_hash{$mkey};
            }
        }
        # p(%this_edge_index_hash);p(@this_nod_index_set);exit;
        ### link index add
        my $before_num=1;
        my $after_num=0;
        while($after_num != $before_num){
            $before_num = $after_num;
            foreach my $skey (keys %this_edge_index_hash){
                my @mutiple_set = @{$this_edge_index_hash{$skey}};
                my $stat_num=0;
                foreach my $ad (@mutiple_set){
                    if (ref($ad) eq 'ARRAY'){
                        my $plus = ${$ad}[1] + 1;
                        my $minus = ${$ad}[0] -1;
                        if (grep {$_ == $plus} @this_nod_index_set ){
                            push @this_nod_index_set,@$ad;
                            $stat_num++;
                        }elsif( grep {$_ == $minus} @this_nod_index_set){
                            push @this_nod_index_set,@$ad;
                            $stat_num++;
                        }elsif (grep {$_ == $ad } @this_nod_index_set){
                            $stat_num++;
                        }
                    }else{
                        my $plus = $ad+1;
                        my $minus = $ad-1;
                        if (grep {$_ == $plus} @this_nod_index_set ){
                            push @this_nod_index_set,$ad;
                            $stat_num++;
                        }elsif( grep {$_ == $minus} @this_nod_index_set){
                            push @this_nod_index_set,$ad;
                            $stat_num++;
                        }elsif (grep {$_ == $ad } @this_nod_index_set){
                            $stat_num++;
                        }
                    }
                }
                if($stat_num > 0){
                    delete $this_edge_index_hash{$skey};
                    $after_num++;
                }
            }
            if (keys %this_edge_index_hash == 0){
                last;
            }
        }
        ### check the barcode type below
        my @this_solution_support;
        # p(%this_edge_index_hash);p(@this_nod_index_set);exit;
        # if ($barcode eq 'barcode_57199'){
            # p(%this_edge_index_hash);p(@this_nod_index_set);
        # }
        if (keys %this_edge_index_hash == 0){
            ### 0-1
            #print OUTB 'This barcode not have mutiple region, next to 0-1 '.$barcode."\n";
            my @fix_nset = sort {$a<=>$b} uniq(@this_nod_index_set);
            my $min_index = min(@fix_nset);
            my $max_index = max(@fix_nset);
            if ($chrtype eq 'cycle'){
                if ($max_index==$max_edge and $min_index==0){
                    my ($r_a,$r_b) = split_two_region(\@fix_nset);
                    unless ($r_a){
                        print OUTB 'Error Don\'t split two region for cycle chrmosome The 0-1 way barcode name '.$barcode."\t".'solution name '.$solution."\n";
                    }else{ 
                    ###region_a
                        my ($a_len,$a_gap) = get_length_and_gap(\@$r_a,\@path_set,$barcode,$solution,'0-1-1');
                        ###region_b
                        my ($b_len,$b_gap) = get_length_and_gap(\@$r_b,\@path_set,$barcode,$solution,'0-1-2');
                        my $this_max_gap = max($a_gap,$b_gap);
                        my $all_length = $a_len+$b_len;
                        @this_solution_support = ($all_length,$this_max_gap);
                    }
                }else{
                    my ($this_barcode_len,$max_gap) = get_length_and_gap(\@fix_nset,\@path_set,$barcode,$solution,'0-2');
                    @this_solution_support = ($this_barcode_len,$max_gap);
                    #print OUTD $barcode."\t".join("\t",@this_solution_support)."\n";
                }
            }else{
                my ($this_barcode_len,$max_gap) = get_length_and_gap(\@fix_nset,\@path_set,$barcode,$solution,'0-2');
                @this_solution_support = ($this_barcode_len,$max_gap);
                #print OUTD $barcode."\t".join("\t",@this_solution_support)."\n";
            }
        }elsif (keys %this_edge_index_hash == 1){
            #print OUTB 'This barcode have one mutiple region, next to 1-1 '.$barcode."\n";
            my @key_set = keys %this_edge_index_hash;
            my @the_index_set = @{$this_edge_index_hash{$key_set[0]}};
            foreach my $adex (@the_index_set){
                my @double_nod_index_set = @this_nod_index_set;
                push @double_nod_index_set,$adex;
                my @fix_nset = sort {$a<=>$b} uniq(@double_nod_index_set);
                my $min_index = min(@fix_nset);
                my $max_index = max(@fix_nset);
                my $b_len;
                my $g_len;
                if ($chrtype eq 'cycle'){
                    if ($max_index == $max_edge and $min_index==0){
                        my ($r_a,$r_b) = split_two_region(\@fix_nset);
                        unless ($r_a){
                            print OUTB 'Error Don\'t split two region for cycle chrmosome The 0-1 way barcode name '.$barcode."\t".'solution name '.$solution."\n";
                        }else{
                            ###region_a
                            my $ta_len;
                            my $ta_gap;
                            if (@$r_a == 1){
                                $ta_len = $region_hash{$path_set[${$r_a}[0]]};
                            }else{
                                ($ta_len,$ta_gap) = get_length_and_gap(\@$r_a,\@path_set,$barcode,$solution,'1-1-1');
                            }
                            ###region_b
                            my $tb_len;
                            my $tb_gap;
                            if (@$r_b == 1){
                                $tb_len = $region_hash{$path_set[${$r_b}[0]]};
                                $tb_gap = 0;
                            }else{
                                ($tb_len,$tb_gap) = get_length_and_gap(\@$r_b,\@path_set,$barcode,$solution,'1-1-2');
                            }
                            $b_len = $ta_len+$tb_len;
                            $g_len = $ta_gap+$tb_gap;
                        }
                    }else{
                        my ($this_barcode_len,$max_gap) = get_length_and_gap(\@fix_nset,\@path_set,$barcode,$solution,'2-2');
                        $b_len = $this_barcode_len;
                        $g_len = $max_gap;
                    }
                }else{
                    my ($this_barcode_len,$max_gap) = get_length_and_gap(\@fix_nset,\@path_set,$barcode,$solution,'1-2');  
                    $b_len = $this_barcode_len;
                    $g_len = $max_gap;
                }
                if (@this_solution_support == 0){
                    @this_solution_support = ($b_len,$g_len);
                }else{
                    if ($b_len < $this_solution_support[0]){
                        @this_solution_support = ($b_len,$g_len);
                    }
                }
            }
        }else{
            ###2-1
            my $duplication_num = scalar(keys %this_edge_index_hash);
            if ($duplication_num > $mbdup){
                #print OUTB 'This barcode more than '.$mbdup.' mutiple region, jump!!!'.$barcode."\n";
                next;
            }
            my @merge_set;
            my %class_set = check_duplication(\%this_edge_index_hash);
            foreach my $keys (keys %class_set){
                push @merge_set,[@{$this_edge_index_hash{$keys}}];
            }
            array_permute(@merge_set);
            foreach my $set (@all_set){
                my @double_nod_index_set = @this_nod_index_set;
                push @double_nod_index_set,@$set;
                my @fix_nset = sort {$a<=>$b} uniq(@double_nod_index_set);
                my $min_index = min(@fix_nset);
                my $max_index = max(@fix_nset);
                my $b_len;
                my $g_len;
                if ($chrtype eq 'cycle'){
                    if ($max_index==$max_edge and $min_index==0){
                        my ($r_a,$r_b) = split_two_region(\@fix_nset,$barcode,$solution,'2-1-1');
                        unless ($r_a){
                            print OUTB 'Error Don\'t split two region for cycle chrmosome The 0-1 way barcode name '.$barcode."\t".'solution name '.$solution."\n";
                        }else{
                            ###region_a
                            my $ta_len;
                            my $ta_gap;
                            if (@$r_a == 1){
                                $ta_len = $region_hash{$path_set[${$r_a}[0]]};
                            }else{
                                ($ta_len,$ta_gap) = get_length_and_gap(\@$r_a,\@path_set,$barcode,$solution,'2-1-1');
                            }
                            ###region_b
                            my $tb_len;
                            my $tb_gap;
                            if (@$r_b == 1){
                                $tb_len = $region_hash{$path_set[${$r_b}[0]]};
                                $tb_gap = 0;
                            }else{
                                ($tb_len,$tb_gap) = get_length_and_gap(\@$r_b,\@path_set,$barcode,$solution,'2-1-2');
                            }
                            $b_len = $ta_len+$tb_len;
                            $g_len = $ta_gap+$tb_gap;
                        }
                    }else{
                        my ($this_barcode_len,$max_gap) = get_length_and_gap(\@fix_nset,\@path_set,$barcode,$solution,'2-2');
                        $b_len = $this_barcode_len;
                        $g_len = $max_gap;
                    }
                }else{
                    #2-2
                    my ($this_barcode_len,$max_gap) = get_length_and_gap(\@fix_nset,\@path_set,$barcode,$solution,'2-2');
                    $b_len = $this_barcode_len;
                    $g_len = $max_gap;
                }
                if (@this_solution_support == 0){
                    @this_solution_support = ($b_len,$g_len);
                }else{
                    if ($b_len < $this_solution_support[0]){
                        @this_solution_support = ($b_len,$g_len);
                    }
                }     
            }
            @all_set=();
        }
        # p(@this_solution_support);p($max_allow_length);exit;
        if (@this_solution_support){
                if ($this_solution_support[0] < $max_allow_length and $this_solution_support[1] < 10 ){
                    my $t_bar_num = ${$bs_set{$barcode}}[2];
                    $bs_stat = $bs_stat + $t_bar_num;
            }   
        }else{
            print OUTB 'Error not calculate the barcode lengt and gap length '."$solution\t$barcode\n";
        }
    }
    print OUTA "$bs_stat\t$solution\n";
}

###########################
sub re_sort_edge{
    my $input = $_[0];
    chomp($input);
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

#######################
sub re_index_edge{
    my ($l_set,$t_set) = @_;
    my %index_set;
    foreach my $ledg (@$t_set){
        my @tdex = indexes{$_ eq $ledg} @$l_set;
        @{$index_set{$ledg}} = @tdex;
    }
    return %index_set;
}
####################################################
sub re_index_nod{
    my ($n_set,$s_set) = @_;
    my @nset;
    foreach my $lnod (@$s_set){
        my @ndex = indexes{$_ eq $lnod} @$n_set;
        if (@ndex){
            push @nset,($ndex[0],$ndex[0]+1);
        }else{
            my $all_num = scalar(@$n_set);
            push @nset,(0,$all_num);
        } 
    }
    return @nset;
}
####################################################
sub stat_length{
    my $apath = $_[0];
    my $the_length = 0;
    foreach my $t_path (@$apath){
        $the_length = $the_length + $region_hash{$t_path};
    }
    return $the_length;
}

#############################################################
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

#############################################################
sub get_path_gap{
    my ($input,$solu_set) = @_;
    my $start_point = shift @$input;
    my @gap_set;
    my $gap_stat=0;
    foreach my $second_point (@$input){
        my $gap = $second_point - $start_point;
        if ($gap > 1){
            my $gap_len=0;
            my @split_set = @$solu_set[int($start_point)+1..int($second_point)-1];
            # my $gap_len = stat_length(\@split_set);
            foreach my $t_edge (@split_set){
                if (grep {$t_edge eq $_} @short_length_edge){
                    $gap_len = $gap_len + 1;
                }else{
                    $gap_len = $gap_len + $region_hash{$t_edge};
                }
            }
            push @gap_set,$gap_len;
            $gap_stat++;
        }
        $start_point = $second_point;
    }
    if ($gap_stat == 0){
        @gap_set = qw/ 1 /;
        return @gap_set;
    }else{
        return @gap_set;
    }
}

##################################################
sub check_duplication{
    my $input = $_[0];
    my %rp_set;
    foreach my $keys (keys %$input){
        my @the_values_set = @{${$input}{$keys}};
        my @double_vals = @the_values_set;
        my @final_set;
        my @the_class_set;
        my $one_point = shift @the_values_set;
        my $leiji=1;
        my $first_index = 0;
        while( @the_values_set !=0){
            my $two_point = shift @the_values_set;
            my $t_gap = $two_point - $one_point;
            if ($t_gap == 1){
                if (@the_values_set == 0){
                    push @final_set,[$one_point,$two_point];
                }
            }else{
                my $decline = $leiji - $first_index;
                if ($decline == 1){
                    push @final_set,$one_point
                }else{
                    push @final_set,[@double_vals[$first_index..$leiji-1]];
                    
                }
                if (@the_values_set == 0){
                    push @final_set,$two_point;
                }
                $first_index = $leiji;
            }
            $one_point = $two_point;
            $leiji++;
        }
        @{$rp_set{$keys}} = @final_set;
    }
    return %rp_set;
}
##############################################
sub get_length_and_gap{
    my ($f_set,$p_set,$ab,$ak,$num)= @_;
    my $min_nset = int(min(@$f_set));
    my $max_nset = int(max(@$f_set));
    if ($min_nset == 0 and $max_nset == 0){
        print OUTB 'Error Length of barcode and gap is empty! '."$ab\t$ak\t$num\n";
    }
    my @split_array = @$p_set[$min_nset..$max_nset];
    my $the_barcode_len = stat_length(\@split_array);
    # print 'sub:the_barcode_len '.$the_barcode_len."\n";
    my @this_gap_set = get_path_gap(\@$f_set,\@$p_set);
    my $m_gap=0;
    if (@this_gap_set){
        $m_gap = max(@this_gap_set);
    }
    return ($the_barcode_len,$m_gap)
}

##############################################
sub split_two_region{
    my $input = $_[0];
    my @copy_set = @$input;
    my $ele_num = scalar(@$input);
    my $first_p = shift @$input;
    my $stat_dex=0;
    my @p_set;
    foreach my $second_p (@$input){
        my $gap = $second_p - $first_p;
        push @p_set,$gap;
    }
    my $m_p = max(@p_set);
    my @g_index_set = indexes{$_ eq $m_p} @p_set;
    my @head = @copy_set[0..int($g_index_set[0])];
    my @tail = @copy_set[int($g_index_set[0]+1)..int($ele_num-1)];
    return (\@head,\@tail);
}

######################################################3
sub array_permute{
    my $array_input = shift @_;
    foreach (@$array_input){
        if (ref($_) eq 'ARRAY'){
            push @arr,@$_;
            my $element = scalar(@$_);
            array_permute(@_) if (@_);
            push @all_set,[@arr] unless (@_);
            for(my $i=0;$i<$element;$i++){
                pop @arr;
            }
        }else{
            push @arr,$_;
            array_permute(@_) if (@_);
            push @all_set,[@arr] unless (@_);
            pop @arr;
        }
    }
}