#!/usr/bin/perl -w

use strict;
use Data::Printer;
use Cwd qw/getcwd/;
use Getopt::Long;
use List::Util qw/max min sum tail/;
use Array::Utils qw(:all);
use List::MoreUtils qw/uniq/;
use Classbarcode qw/class_barcode same_class_barcode/;
use POSIX qw/floor/;

print 'Start time:'.localtime()."\n";
my $usage=<<USAGE;

    v16 pangwending 2022.06.01
    perl $0 <barcode_stat> <nod> <edg> <total> [-option]

        -n              :Sample name
        -o              :outdir
        -chr_type       :[liner | cycle] the chrmosome type 
        -fix            :[undef | others]Pleas input the new duplication type! [example: 14_15,T,T;8_9,T,F (edge,tandem,untandem)]
        -b              :the link-edge barcode use max [N+b] [defacult:b=1]
        -fp             :[auto | others number] the number of the link-path barcode support. [defacult:auto]
        -i              :if produce the complete genome result(info file) and [default:F]
        -len            :the cut length of the path to find barcode [defacutl:5000] 
        -start_nod      :[undef | others] the start link nod [defacult:auto-find uniq edge] (manual setup: 1/2,3/4,5/6); 
        -way            :[only |all ] only produce scaffold or continue to filter the complete result [only or all]
 
USAGE


my ($barcode_stat,$nod,$edge,$region) = @ARGV;
my ($name,$outdir,$chr_type,$start_nod,$way,$ainfo,$start_length,$repair,$bar_edge_num,$f_path);

GetOptions(
    "n:s" => \$name,
    "i:s" => \$ainfo,
    "len:s" => \$start_length,
    "fix:s" => \$repair,
    "b:s" => \$bar_edge_num,
    "fp:s" => \$f_path,
    "chr_type:s" => \$chr_type,
    "start_nod:s" => \$start_nod,
    "way:s" =>\$way,
    "o:s" => \$outdir
);

$ainfo ||='F';
$name ||='test';
$way ||='only';
$start_length ||='5000';
$bar_edge_num ||='1';
$f_path ||='auto';
$outdir ||=getcwd;
$start_nod ||='undef';
$repair ||='undef';

die $usage if (!$barcode_stat||!$nod||!$edge||!$chr_type);

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
my $origin_edge;
my @uniq_edge_set;
my $j=1;
open INA,$edge or die $!;
while(<INA>){
    chomp;
    next if (/^\#/);
    my ($edge,$cn) = (split,/\s+/,$_)[0,1];
    next if ($cn == 0);
    $edge_cn_hash{$edge} = $cn;
    $origin_edge = $edge if ($j==1);
    push @uniq_edge_set,$edge if ($cn == 1);
    my ($l,$r) = split/_/,$edge;
    $edge_index_hash{$l} = ($edge);
    my $b_edge = $r.'_'.$l;
    $edge_index_hash{$r} = ($b_edge);
    $j++;
}
close INA;

##### nod set file information ###################
open INB,$nod or die $!;
my %nod_all;
my %nod_hash_set;
while(<INB>){
    chomp;
    my ($a_n,$num) = split/\s+/,$_;
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

##### stat the useful barcode number ################
my $useful_barcode_num = int(system("cat $barcode_stat |wc -l"));
my $median = 200;
if ($f_path eq 'auto'){
    if ($useful_barcode_num < 200){
        $f_path = 3;
    }else{
        $f_path = 3;
        for (my $i=200;$i<$useful_barcode_num;$i=$i+200){
            if ($useful_barcode_num - 200 < 200){
                $f_path = $f_path + 3;
                last; 
            }else{
                $f_path = $f_path + 3;
            }
        }
    }
}

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

###############1. get the start nod ###################
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
        my $rever_this_path = link_reverse($this_path);
        push @c_set,($this_path,$rever_this_path);
    }else{
        print '!left and right_edge is not exists! This is not uniq edge!'."\n";
    }
}
# p(@c_set);exit;
############2.complete the nod path ###############
my %path_set;
if (@c_set){
    foreach my $line (@c_set){
        $path_set{$line} = '1';
    }
}else{
    if ($start_nod ne 'undef'){
        my @alon_nod = split/\,/,$start_nod;
        foreach my $anod (@alon_nod){
            my @index = split/\//,$anod;
            my $left_link = judge_odd_even($index[0],'l');    
            my $right_link = judge_odd_even($index[1],'r');
            my $path = $left_link.'/'.$right_link;
            my ($l,$r) = split/_/,$left_link;
            $path_set{$path} = '1';
        }
    }else{
        print '!Not exists uniq edge Please order the nod! You can manual setup the start_nod set!'."\n";
        exit;
    }
}

#############3. use barcode information to produce scaffold ###############
my @final_set;
my %sub_path = %path_set;
my $this_genome_length = sum(values %edge_cn_hash);
my $k=0;
###########3.make link by the special barcode ############
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
        if ($this_uniq_edge_number == 2){
            @end_two_edge = @this_edge_set[-2,-1];
            @last_nod = ($this_nod_set[-1]);
        }else{
            my @to_stat_length = @this_edge_set[-2,-1];
            my $two_edge_length = stat_length(\@to_stat_length,\%region_hash);
            if ($two_edge_length < $start_length){
                my @cut_set = cut_order_length(\@this_edge_set,3);
                if (@cut_set){
                    my $nodnum = pop @cut_set;
                    @end_two_edge = @cut_set;
                    @last_nod = tail $nodnum-1,@this_nod_set;
                }else{
                    @end_two_edge = @this_edge_set[-2,-1];
                    @last_nod = ($this_nod_set[-1]);
                }
            }else{
                @end_two_edge = @this_edge_set[-2,-1];
                @last_nod = ($this_nod_set[-1]);
            }
        }
        my $last_point = (split/_/,$first)[-1]; 
        if (@this_edge_set == $this_genome_length){
            push @final_set,$first;
            delete $sub_path{$first};
            print '#This path genome length is enough! delete and output this path:'.$first."\n";
            next;
        }
        ###check left link way ###############
        my $this_new_path = check_path_nod($first,$last_point,\%nod_all);
        #p($this_new_path);exit;
        if ($this_new_path ne 'undef'){
            if ($this_new_path eq 'BK'){
                push @final_set,$first;
                delete $sub_path{$first};
                print '#This point is the head point and tail point. End the link!'."\n";
                next;
            }
            if (control_the_cnv($this_new_path,\%edge_cn_hash) eq 'T'){
                $sub_path{$this_new_path} = '1'; 
                delete $sub_path{$first};
                print '#This path have a only link-nod! delete this path: '.$first."\n";
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
                    my %new_barcode_hash = class_barcode($barcode_stat,\@end_two_edge,\@last_nod,$bar_edge_num);
                     
                    my %this_bs_stat = untandem_link($first,\%new_barcode_hash,\@last_nod,$last_point,\%edge_cn_hash,$copy_point);
                    my @the_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                    print '#This sub_path location is untandeme!'."\n";
                    if (exists $the_link_set[0]){
                        my @uniq_set = uniq (@the_link_set);
                        if (@uniq_set == 1){
                            $sub_path{$uniq_set[0]} = '1';
                            print '#We get the possible link path: '.$uniq_set[0]."\n";
                        }else{
                            push @final_set,$first;
                            print '#Mutiple link result. link failed!!!'."\n";
                            print join(",",@uniq_set)."\n";
                            print '#This barcode set:'."\n";
                            print_hash(\%new_barcode_hash);
                            print_hash(\%this_bs_stat);
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
                                print '#We get this cycle final link path: '.$uniq_set[0]."\n\n";
                            }else{
                                push @final_set,$first;
                                print '#Mutiple link result. link failed!!!'."\n";
                                print join(",",@uniq_set)."\n";
                                print '#This barcode set:'."\n";
                                print_hash(\%this_barcode_hash);
                                print_hash(\%this_bs_stat);
                            }
                        }else{
                            push @final_set,$first;
                            print '#No result! link failed!!!'."\n";
                            print '#This barcode set:'."\n";
                            print_hash(\%this_barcode_hash);
                        }
                    }
                ### 2.1.3 
                }elsif($tandem_num > 0 and $untandem_num == 0){
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
                            print '#We get this cycle final path:'.$uniq_set[0]."\n\n";
                        }else{
                            print '#Mutiple link result. link failed!!!'."\n";
                            print join(",",@uniq_set)."\n";
                            print '#This barcode set:'."\n";
                            print_hash(\%this_barcode_hash);
                            print_hash(\%this_bs_stat);
                            push @final_set,$new_path;
                        }
                    }else{
                        print '#No result. link failed!!!'."\n";
                        print '#This barcode set:'."\n";
                        print_hash(\%this_barcode_hash);
                        push @final_set,$new_path;
                    }  
                }else{
                    print '!exists the tandem and untandem barcode? This barcode can\'t know this lcation infroamtion!'."\n\n";
                    print '#This path is '.$first."\n";
                    print '#This barcode set:'."\n";
                    print_hash(\%this_barcode_hash);
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
                print '#We get the complete tandem duplication path:'.$link_result."\n";
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
                        print '#We get this cycle final path:'.$uniq_set[0]."\n\n";
                    }else{
                        push @final_set,$link_result;
                        print '#Mutiple link result. link failed!!!'."\n";
                        print join(",",@uniq_set)."\n";
                        print '#This barcode set:'."\n";
                        print_hash(\%inthis_t_bs);
                        print_hash(\%this_bs_stat);
                    }
                }else{
                    push @final_set,$link_result;
                    print '#No link result. link failed!!!'."\n";
                    print '#This barcode set:'."\n";
                    print_hash(\%inthis_t_bs);
                }
            ### 2.3
            }elsif($dupli_type_set[0] eq 'F' and $dupli_type_set[1] eq 'T'){
                print '#This sub_path next link is the F-T'."\n";
                my %this_barcode_hash = class_barcode($barcode_stat,\@end_two_edge,\@last_nod,$bar_edge_num);
                my %this_bs_stat = untandem_link($first,\%this_barcode_hash,\@last_nod,$last_point,\%edge_cn_hash);
                my @the_link_set = barcode_support_num(\%this_bs_stat,$f_path);
                if (exists $the_link_set[0]){
                    my @uniq_set = uniq (@the_link_set);
                    if (@uniq_set == 1){
                        $sub_path{$uniq_set[0]} = '1';
                        print '#We get this cycle final path:'.$uniq_set[0]."\n\n";
                    }else{
                        push @final_set,$first;
                        print '#Mutiple link result. link failed!!!'."\n";
                        print join(",",@uniq_set)."\n";
                        print '#This barcode set:'."\n";
                        print_hash(\%this_barcode_hash);
                        print_hash(\%this_bs_stat);
                    }
                }else{
                    push @final_set,$first;
                    print '#No link result. link failed!!!'."\n";
                    print '#This barcode set:'."\n";
                    print_hash(\%this_barcode_hash);
                }
            }else{
                print '#It\'s impossible for the duplication type!'."\n";
            }
        ### 2. the alon edge link ##########################
        }else{
            print '#This link only one-nod:'."\n";
            my @cnod = @{$nod_all{$last_point}};
            my $rdex = (split/\//,$cnod[0])[1];
            my $next_edge = $edge_index_hash{$rdex};
            my $new_path = $first.'/'.$next_edge;
            $sub_path{$new_path} = '1';
            print '#We get new path:'.$new_path."\n\n";
        }
        delete $sub_path{$first} if (exists $sub_path{$first});
    }
    $k++;
    last if (keys %sub_path == 0);
}

###### rm the duplication path and merge the common path #########
print "!!!Finally $k cycle in this!\n\n";

my @output_set;
foreach my $this_solution (@final_set){
    my $this_edge_number = split/\//,$this_solution;
    if ($this_edge_number == $this_genome_length){
        my ($edge_confict,$nod_confict) = control_confict($this_solution,\%edge_cn_hash,\%nod_hash_set);
        if ($edge_confict+$nod_confict == 0){
            push @output_set,$this_solution;
        }
    }
}

if ($ainfo eq 'T'){
    open OUTC,">$outdir/$name.info" or die $!;
}

open OUTB,">$outdir/$name.result" or die $!;
if (@output_set){
    print OUTB join("\n",@output_set)."\n";
    print '#exists complete path! This link is end.'."\n";
}else{
    print '#Before the rmdup scaffold set: '.join("\t",@final_set)."\n";
    my @rm_final_set;
    if ($chr_type eq 'cycle'){
        my $min_dex = min(keys %nod_all);
        my @cut_head_tail = cut_cycle_zero(\@final_set,$min_dex);
        @rm_final_set = rmduplicate(\@final_set);
    }else{
        @rm_final_set = rmduplicate(\@final_set);
    }
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
    print '#Next start to get the all genome link path set! loading........'."\n";
    my @result = right_direct_link($origin_edge,\%edge_index_hash,\%nod_all,\%edge_cn_hash,\@rm_final_set,$ainfo,\%nod_hash_set);
    # p(@result);exit;
    if ($ainfo eq 'T'){
        close OUTC;
    }
    if (@result){
        print '#Output the final screen set!'."\n";
        print OUTB join("\n",@result)."\n";
        close OUTB;
    }else{
        print '!It\'s so pity! No result output. You should check the barcode file repair the parameter or use the info set!.'."\n";
        close OUTB;
    }
}

print '****************************This work is finished!!!*************************************'."\n";
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
    }else{
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
                        push @{$bs_stat{$link_result}},'1';
                    }
                }
            }
        }
    }
    return %bs_stat;
}
#############sub right_direct_link ################ 
sub right_direct_link{
    my ($origin,$edge_set,$nod_set,$edge_hash,$local_set,$oinfo,$anodset)=@_;
    my $genome_length = sum(values %$edge_hash);
    my @final_solution;
    my @list = ($origin);
    my $all_path_num = 0;
    my $filter_path_num = 0;
    while(1){
        my $edge1 = shift @list;
        my $this_solution_length = split/\//,$edge1;
        if ($this_solution_length == $genome_length){
            my ($edge_c,$nod_c) = control_confict($edge1,\%$edge_hash,\%$anodset);
            if ($edge_c+$nod_c == 0){
                $all_path_num++;
                if ($oinfo eq 'T'){
                    print OUTC $edge1."\n";
                }
                my @path_corret;
                foreach my $this_local (@$local_set){
                    if ($edge1 =~ /$this_local/){
                        push @path_corret,'T';
                    }else{
                        my $rever_edge = link_reverse($this_local);
                        if ($edge1 =~ /$rever_edge/){
                            push @path_corret,'T';
                        }else{
                            push @path_corret,'F';
                        }
                    }
                }
                unless (grep {$_ eq 'F'} @path_corret){
                    push @final_solution,$edge1;
                    # p($edge1);
                    $filter_path_num++;
                }
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
    print '*All link path number:'.$all_path_num."\t".'Filter link path number:'.$filter_path_num."\n";
    return @final_solution;
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
    foreach my $re_ed (@uniq){
        my $num = grep {$_ eq $re_ed} @set;
        if ($num <= ${$this_hash}{$re_ed}){
            push @stat,'1'; 
        }else{
            push @stat,'0';
        }
    }
    if (grep {$_ eq '0'} @stat){
        return 'F';
    }else{
        return 'T';
    }
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
    my $input = $_[0];
    if (%$input){
        foreach my $keys (keys %$input){
            my @this_set = @{${$input}{$keys}};
            print "$keys\t".join("\t",@this_set)."\n";
        }
    }else{
        print 'Empty!'."\n";
    }
    print "\n";
}

#################sub control_confict #######################
sub control_confict{
    my ($edg,$edge,$nod_set) = @_;
    my @nodset = split/_/,$edg;
    my $first_poin = $nodset[0]; 
    shift @nodset;
    my $last_poin = $nodset[-1];
    pop @nodset;
    my @new_edge_set = re_sort_edge($edg);
    my @erro_edge;
    ##stat the edge confitct
    my %double_edge = %$edge;
    my @uniq_arr = uniq(@new_edge_set);
    foreach my $de2 (@uniq_arr){
        my $num = grep {$_ eq $de2} @new_edge_set;
        if ($num == ${$edge}{$de2}){
            push @erro_edge,'0';
        }elsif($num < ${$edge}{$de2}){
            my $err_num = ${$edge}{$de2} - $num;
            push @erro_edge,$err_num;
        }else{
            my $err_num = $num - ${$edge}{$de2};
            push @erro_edge,$err_num;
        }
        delete $double_edge{$de2} if (exists $double_edge{$de2});
    }
    my @edgeleft = values %double_edge;
    push @erro_edge,@edgeleft;
    my $edge_confict = sum(@erro_edge);
    ##stat the nod confict ##
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
    return ($edge_confict,$nod_confict);
}

################################################################
sub cut_order_length{
    my ($input,$wtnumber) = @_;
    my $edge_num = scalar(@$input);
    my $uniq_path_number = scalar(uniq(@$input)); 
    my @final_set;
    if ($wtnumber > $uniq_path_number){
        return @final_set;
        last;
    }
    my $j;
    for ($j=$wtnumber;$j<=$edge_num;$j++){
        my @cut_set = tail $j,@$input;
        my @uniq_cut_set = uniq(@cut_set);
        if (@uniq_cut_set == $wtnumber){
            @final_set = @uniq_cut_set;
            last;
        }
    }
    push @final_set,$j;
    return @final_set;
}

##################################################
sub cut_cycle_zero{
    my ($input,$min) = @_;
    my @output_set;
    foreach my $ap (@$input){
        chomp($ap);
        if ($ap =~ /\/$min\_/){
            my ($left,$right) = (split/\/$min\_/,$ap,2);
            my $b_right = $min.'_'.$right;
            push @output_set,($left,$b_right);
        }elsif($ap =~ /\_$min\//){
            my ($left,$right) = (split/\_$min\//,$ap,2);
            my $b_left = $left.'_'.$min;
            push @output_set,($b_left,$right);
        }else{
            push @output_set,$ap;
        }
    }
    return @output_set;
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
        my $this_set_num = scalar(@{${$input}{$keys}});
        if ($this_set_num < $filter){
            next;
        }else{
            push @return_set,$keys;
        }
    }
    return @return_set;
}

