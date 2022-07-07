#!/usr/bin/perl -w

use strict;
use Data::Printer;
use POSIX qw/floor/;
use Cwd qw/getcwd/;
use Getopt::Long;

my $map = $ARGV[0];
my ($chr_id,$name,$samtools,$outdir);

my $usage=<<USAGE;

Description: v1.0 2022.04.09 pwd 
             This script is to spilt the reads which unmapping and include loxpsym!

Usage: perl $0 <map> [-option]

        -n          :the file name
        -s          :the path of samtools 
        -chr        :the chrmosome name [must]
        -o          :the output directory

USAGE

GetOptions(
    "n:s" =>\$name,
    "chr:s" =>\$chr_id,
    "s:s" =>\$samtools,
    "o:s" =>\$outdir
);

die $usage if (!$map||!$chr_id);

my $loxpsym = 'ATAACTTCGTATAATGTACATTATACGAAGTTAT';
$samtools ||='/hwfssz1/ST_BIOCHEM/P18Z10200N0255/PMO/F14ZQSYJSY1726/pangwending/some_tools/samtools-1.13/bin/samtools';
$name ||= 'test';
$outdir ||=getcwd;

if ($map =~ /\.sam$/){
    open INC,$map or die $!;
}else{
    open INC,"$samtools view $map $chr_id |" or die $!;
}

my %reads_set;
my $unmap_loxp_reads=0;
open OUTA,">$outdir/$name.split.fq1" or die $!;
open OUTB,">$outdir/$name.split.fq2" or die $!;
open OUTC,">$outdir/$name.loxp.fq" or die $!;
while(<INC>){
    chomp;
    next if (/^\@/);
    my ($id,$flag,$chrid,$start,$map_quality,$align,$seqence,$qual)=(split/\s+/,$_)[0,1,2,3,4,5,9,10];
    my @flag_index = tentotwo($flag);
    if ($flag_index[8]){
        print "$_\n";
        next if ($flag_index[8] == '1');
    }
    unless ($flag_index[2] == '1' or $align eq '*'){
        if ($map_quality >= 30 and $align eq '100M'){
            next;
        }else{
            if ($seqence =~ /$loxpsym/gi){
                $unmap_loxp_reads++;
                my ($split_set1,$split_set2) = split_seqence($id,$seqence,$qual,$loxpsym);
                print OUTA join("\n",@$split_set1)."\n";
                print OUTB join("\n",@$split_set2)."\n";
                print OUTC '@'.$id."\n".$seqence."\n".'+'."\n".$qual."\n";
            }
        }
    }else{
        if ($seqence =~ /$loxpsym/gi){
            $unmap_loxp_reads++;
            my ($split_set1,$split_set2) = split_seqence($id,$seqence,$qual,$loxpsym);
            print OUTA join("\n",@$split_set1)."\n";
            print OUTB join("\n",@$split_set2)."\n";
            print OUTC '@'.$id."\n".$seqence."\n".'+'."\n".$qual."\n";
        }
    }

}
close INC;
close OUTA;
close OUTB;
close OUTC;

######################################################
sub split_seqence{
    my ($my_id,$my_seq,$my_quali,$loxp) = @_;
    my $index = index($my_seq,$loxp);
    my $reads1 = substr($my_seq,0,$index+34);
    my $quali1 = substr($my_quali,0,$index+34);
    my $true_id1 = '@'.$my_id.'/1';

    my $reads2 = substr($my_seq,$index);
    my $quali2 = substr($my_quali,$index);
    my $true_id2 = '@'.$my_id.'/2';
    my @set1;
    my @set2;
    my $len1 = length($reads1)-34 +1;
    my $len2 = length($reads2) -34 + 1;
    if ($len1 >=15 and $len2 >=15){
        push @set1,($true_id1,$reads1,'+',$quali1);
        push @set2,($true_id2,$reads2,'+',$quali2);
    }
    return (\@set1,\@set2);
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


