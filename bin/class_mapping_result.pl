#!/usr/bin/perl -w

use strict;
use Cwd qw/getcwd/;
use List::MoreUtils qw/uniq/;
use Getopt::Long;
use POSIX qw/floor/;

my $usage=<<USAGE;

Description: This script is to stat the alignment result and prepare the data for after use!
             2022.05.11 v1.0  pangwending
Usage: perl $0 <map> [-option]

            -n      :Sample name
            -chr    :chr_id [default:IXR_BACseq]
            -q      :alignment quality (0~60) [default:30]
            -s      :the path of samtools
            -o      :outdir
            -h      :help

Example:   perl $0 -n JS599 -chr IXR_BACseq -m1 JS599.sort.bam -q 30 -class bam -o ./

USAGE

my $map1 =$ARGV[0];
my ($chr,$name,$quality,$samtools,$outdir,$help);
GetOptions(
    "n:s" => \$name,
    "chr:s" => \$chr,
    "q:s"=> \$quality,
    "s:s"=> \$samtools,
    "o:s" => \$outdir,
    "h:s" => \$help
);

$outdir ||=getcwd;
$quality ||='30';
$name ||='test';
$samtools ||= '/hwfssz1/ST_BIOCHEM/P18Z10200N0255/PMO/F14ZQSYJSY1726/pangwending/some_tools/samtools-1.13/bin/samtools';

die $usage if (!$map1||!$name||!$chr||$help);

print "Start time:".localtime()."\n";
print "#stats the alignment result will start,it will make a long time, loading........\n";
my $this_start = time();
if ($map1 =~ /\.sam$/){
    open INA,$map1 || die $!;
}else{
    open INA,"$samtools view $map1 $chr |" || die $!;
}
open OUTA,">$outdir/$name.split.fq" or die $!;
open OUTD,">$outdir/$name.loxp.fq" or die $!;
my $loxpsym = 'ATAACTTCGTATAATGTACATTATACGAAGTTAT';
my $index1=0;
my $loxp_reads_num=0;
my $reads_num=0;
my $split_reads_num=0;
my %syn_pair_hash_set;
my %barcode_num;

while(<INA>){
    chomp;
    next if (/^\@/);
    my ($id,$flag,$chrid,$start,$map_quality,$align,$seqence,$qual)=(split/\s+/,$_)[0,1,2,3,4,5,9,10];
    my $barcode;
    if ($id =~ /\#(.+)/){
        $barcode = $1;
        push @{$barcode_num{$barcode}},'1';
    }
    $reads_num++;
    my @index =  tentotwo($flag);
    if ($index[8]){
        if ($index[8] == '1'){
            if ($seqence =~ /$loxpsym/gi){
                print OUTD '@'.$id."\n".$seqence."\n".'+'."\n".$qual."\n";
                $loxp_reads_num++;
                my @split_set = split_seqence($id,$seqence,$qual,$loxpsym);
                print OUTA join("\n",@split_set)."\n";
                $split_reads_num++;
            }
            next;
        }
    }
    unless ($index[2] == '1'){
        my $end = $start + 100 -1;
        my $dire;
        if ($index[4] == '1'){
            $dire = '-';
        }else{
            $dire = '+';
        }
        if ($map_quality >= $quality and $align eq '100M'){
            push @{$syn_pair_hash_set{$id}},($chrid,$start,$end,$dire);
            $index1++;
        }else{
            if ($seqence =~ /$loxpsym/gi){
                print OUTD '@'.$id."\n".$seqence."\n".'+'."\n".$qual."\n";
                $loxp_reads_num++;
                my @split_set = split_seqence($id,$seqence,$qual,$loxpsym);
                print OUTA join("\n",@split_set)."\n";
                $split_reads_num++;
            }
        }
    }else{
        if ($seqence =~ /$loxpsym/gi){
            print OUTD '@'.$id."\n".$seqence."\n".'+'."\n".$qual."\n";
            $loxp_reads_num++;
            my @split_set = split_seqence($id,$seqence,$qual,$loxpsym);
            print OUTA join("\n",@split_set)."\n";
            $split_reads_num++;
        }
    }        
}
close OUTA;
close OUTD;

open OUTC,"| gzip >$outdir/$name.synmap.gz" or die $!;
foreach my $keys1 (keys %syn_pair_hash_set){
    my @set1 = @{$syn_pair_hash_set{$keys1}};
    next if (@set1 < 8);
    print OUTC "$keys1\t".join("\t",@set1)."\n";
}

close OUTC;
my $barcode_number = keys %barcode_num;
print '#The numbers of reads: '.$reads_num."\n";
print '#The barcode number: '.$barcode_number."\n";
print '#The numbers of syn chr reads: '.$index1."\n";
print '#The number of loxpsym reads: '.$loxp_reads_num."\n";
print '#The number of split-reads num: '.$split_reads_num."\n";

print "#This work is finshed!\n";
print "End:".localtime()."\n";
my $end_time = time();
my $use_time = ($end_time-$this_start)/60;
print "$use_time\n";
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
    my $len1 = length($reads1)-34 + 1;
    my $len2 = length($reads2)-34 + 1;
    my @set;
    if ($len1 >=15 and $len2 >=15){
        push @set,($true_id1,$reads1,'+',$quali1);
        push @set,($true_id2,$reads2,'+',$quali2);
    }
    return @set;
}