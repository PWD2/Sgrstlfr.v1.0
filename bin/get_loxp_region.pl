#!/usr/bin/perl -w

use strict;
#use Data::Printer;

my ($fasta,$chrid,,$name,$outdir,$start_p) = @ARGV;

my $usage=<<USAGE;

Description: v1.0 2022 08 01 pangwending wangyun

             This script is get the loxp basic information in the syn-chr fasta file !
             output the file inculde [test.list.loxpreg] and [test.lxpftr.lst]

Usage: perl $0 <fasta> <chrid/list> <name> <outdir> <start_p>

USAGE

die $usage if (!$fasta);

my $loxpsym = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";

open OUTA,">$outdir/$name.list.loxpreg" or die $!;
open OUTB,">$outdir/$name.lxpftr.lst" or die $!;

my @chr_list;
if ($chrid =~ /\,/){
    @chr_list = split/\,/,$chrid;
}else{
    push @chr_list,$chrid;
}
my $loxplen = length($loxpsym);
my %seqence = read_fasta($fasta,\@chr_list);


$start_p ||=1;
my $loxp_index = 1;
foreach my $keys (keys %seqence){
    my $index = 0;
    my $start = 0;
    my $end = 1;
    my @loxp_reg;
    
    my $this_seq = $seqence{$keys};
    my $seq_len = length($this_seq);
    while($index != -1){
        $index = index($this_seq,$loxpsym,$start);
        if ($index != -1){
            $start = $index + 17 ;
            my $l_start = $index + 1;
            my $l_end = $l_start+$loxplen -1;
            print OUTA "$keys\tloxp\tloxpreg\t$end\t$start\t.\t.\t.\t$loxp_index\n";
            print OUTB "FTR\tloxp\t$keys\t$l_start\t$l_end\n";
            $end = $start + 1;
            $loxp_index++;
        }else{
            $start = $start+1;
            print OUTA "$keys\tloxp\tloxpreg\t$start\t$seq_len\t.\t.\t.\t$loxp_index\n";
        }
    }
    if ($loxp_index % 2 == 0){
        $loxp_index = $loxp_index + 3;
    }else{
        $loxp_index = $loxp_index + 2;
    }
    
}
open OUTC,"$outdir/$name.ref.encod" or die $!;
my @index_set;
foreach my $i (1..$loxp_index-1){
    push @index_set,$i;
}
print OUTC "$name\t".join(",",@index_set)."\n";
close OUTA;
close OUTB;
print '!This work is finished!'."\n";

########## read_fasta ########################
sub read_fasta{
    my ($input,$chr)=@_;
    $/=">";
    open IND,$input or die $!;
    <IND>;
    my %hash_set;
    while(<IND>){
        chomp;
        my ($i_d,$s_eq)=(split/\n/,$_,2)[0,1];
        if (grep {$i_d eq $_ } @$chr){
            $s_eq =~ s/\r/\n/g;
            $s_eq =~ s/\n//g;
            my $syn_seq = $s_eq;
            $hash_set{$i_d} = $syn_seq;
        }
    }
    $/="\n";
    if(keys %hash_set == 0){
        print '! Not find the appoint chrmosome name seqence. Please corret the fasta file and chrmosome name!'."\n";
        exit;
    }
    close IND;
    return %hash_set;
}

