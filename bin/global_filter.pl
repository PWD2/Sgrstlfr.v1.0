#!/usr/bin/perl -w

use strict;
use Data::Printer;
use List::Util qw/max/;
use List::MoreUtils qw/uniq firstidx/;
use Array::Utils qw(:all);
use Getopt::Long;
use Cwd qw/getcwd/;

print 'Start time:'.localtime()."\n";
my $usage=<<USAGE;

    Description: This script is filter the mutiple assembly result by barcode information! 
    Usage: perl $0 <mutiple-result> <info> <barcode> [-option]

        -n              :file name
        -o              :outdir
        -min            :the minnum of the barcode include edge [default:3]
        -max            :the maxnum of the barcode include edge [default:7]
        -wtnod          :the number of diff nod for stat barcode support [default: 2]
        -way            :the filter way for the barcode support num( max or plus) [default:max]        

USAGE

my ($mutiple,$info,$barcode)=@ARGV;
my ($minedge,$maxedge,$wtnod,$way,$name,$outdir);

GetOptions(
    "n:s" => \$name,
    "o:s" => \$outdir,
    "min:s" => \$minedge,
    "max:s" => \$maxedge,
    "wtnod:s" => \$wtnod,
    "way:s" => \$way
);

die $usage if (@ARGV != 3);
$name ||='test';
$outdir ||=getcwd;
$minedge ||='3';
$maxedge ||='7';
$wtnod ||='2';
$way ||='max';

unless (-s $mutiple){
    print '!!!The mutiple result is Empty!'."\n";
    print 'It will use the info file to filter !'."\n";
    unless (-s $info){
        print 'But the info file is not exists !!! End !'."\n";
        exit;
    }
}
my %barcode_hash;
open INA,$barcode or die $!;
while(<INA>){
    chomp;
    my @this_set = split/\s+/,$_;
    if ($this_set[1] <= $maxedge and $minedge <= $this_set[1]){
        my @s_nod = split/\,/,$this_set[4];
        push @{$barcode_hash{$this_set[0]}},(@s_nod);
    }
}
close INA;

my %path_nod_set;
my $nod_number;
if (-s $mutiple){
    open INB,$mutiple or die $!;
}else{
    open INB,$info or die $!;
}

while(<INB>){
    chomp;
    my @nod_set = re_sort_nod($_);
    $nod_number = scalar(@nod_set);
    @{$path_nod_set{$_}} = @nod_set;
}
close INB;

my $i=0;
open OUTA,">$outdir/$name.final" or die $!;
while($i<$nod_number){
    my @line_nod = get_appoint_one_nod(\%path_nod_set,$i);
    my @uniq_nod = uniq(@line_nod);
    if (@uniq_nod != 1){
        my %line_stat;
        print '#This diff nod set: '.join("\t",@uniq_nod)."\n";
        foreach my $keys1 (keys %path_nod_set){
            my @this_nod_set = @{$path_nod_set{$keys1}};
            my $bar_support = get_dff_place_nod(\@uniq_nod,\@this_nod_set,$wtnod,\%barcode_hash);
            $line_stat{$keys1} = $bar_support;
            print 'Path: '.$keys1."\t".'BS: '.$bar_support."\n";
        }
        %path_nod_set = delete_mutiple(\%line_stat,\%path_nod_set,$way);
    }
    if (keys %path_nod_set == 0){
        print '!!!The all solution is delete! The mutiple maybe wrong!!!'."\n";
        last;
    }
    if (keys %path_nod_set == 1){
        my @output = keys %path_nod_set;
        print '#Left the alon path! Output the end solution!'."\n";
        print OUTA join("\n",@output)."\n";
        last;
    }
    $i++;
    if ($i == $nod_number){
        print '#The nod ergodic is end! Output the end solution!'."\n";
        my @output = keys %path_nod_set;
        print OUTA join("\n",@output)."\n";
    }
}
close OUTA;
print '****************************This work is finished!!!*************************************'."\n";
print 'End time:'.localtime()."\n";
################################################
sub get_appoint_one_nod{
    my ($list,$dex) = @_;
    my @line;
    my $line_num = 1;
    foreach my $keys (keys %$list){
        my @a_set = @{${$list}{$keys}};
        push @line,$a_set[$dex];
    }
    return @line;
}

###############################################
sub get_dff_place_nod{
    my ($bk_set,$anod_set,$wana,$barcode_infor) = @_;
    my $support=0;
    my @the_bk_set;
    my $a_nod_num = scalar(@$anod_set);
    foreach my $bk (@$bk_set){
        my $this_index = firstidx{$_ eq $bk} @$anod_set;
        my $left_dex = $this_index - $wana;
        my $right_dex = $this_index + $wana;
        if ($left_dex < 0){
            $left_dex = 0;
        }
        if ($right_dex >= $a_nod_num){
            $right_dex = $a_nod_num-1;
        }
        my @aset = @$anod_set[$left_dex..$right_dex];
        my $this_num = find_barcode(\@aset,\%$barcode_infor);
        $support = $support + $this_num;
    }
    return $support;
}

###############################################
sub find_barcode{
    my ($path_bk_set,$bar_infor) = @_;
    my $stat_num = 0;
    my $nod_num = scalar(@$path_bk_set);
    # print 'filter set:'.join("\t",@$path_bk_set)."\n";
    foreach my $keys2 (keys %$bar_infor){
        my @bar_nod_set = @{${$bar_infor}{$keys2}};
        if (intersect(@bar_nod_set,@$path_bk_set) == $nod_num){
            # print 'barcode set:'.join("\t",@bar_nod_set)."\n\n";
            $stat_num++;
        }
    }
    return $stat_num;
}

##################################################################
sub re_sort_nod{
    my $input =$_[0];
    my @nod_set = split/_/,$input;
    shift @nod_set;
    pop @nod_set;
    my @new_set;
    foreach my $this_nod (@nod_set){
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

##################################################################
sub delete_mutiple{
    my ($input,$all_path,$get_way) = @_;
    my %new_set = %$all_path;
    if ($get_way eq 'max'){
        my @this_valu = values %$input;
        my $max = max(@this_valu);
        print '#Remove the BS not equal max!'."\n";
        foreach my $this_key (keys %$input){
            my $this_solu_num = ${$input}{$this_key};
            unless ($this_solu_num == $max){
                delete $new_set{$this_key};
                print 'Delete solution: '.$this_key."\n";
            }
        }
    }elsif($get_way eq 'plus'){
        print '#Remove the BS equal Zero!'."\n";
        foreach my $this_key (keys %$input){
            my $this_solu_num = ${$input}{$this_key};
            if ($this_solu_num == 0){
                delete $new_set{$this_key};
                print 'Delete solution: '.$this_key."\n";
            }
        }
    }else{
        print 'Pleas input the corret [way] paramater!'."\n";
        exit;
    }
    print "\n";
    return %new_set;
}

