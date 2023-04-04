#!/usr/bin/perl -w

use strict;
use List::Util qw/pairs/;
use Getopt::Long;
#use Data::Printer;
use FindBin qw($Bin);
use Cwd;

my $usage=<<USAGE;

    Description: run more Sgrstlfr from alist which contain the sample_name,fastq_directory,fastq1,fastq2
    perl $0 <list> <fasta> <refcoverage> [-option]

        <list>          : sample_name,fastq_directory,fastq1,fastq2
        --chrtype       : [linear | cycle ]the chrmosome type (linear or cycle) [default: cycle] 
        --chrid         : the chrid of wanna to restructure [default:IXR_BACseq]
        --t             : the threads about this script run [default: 4]
        --p             : if produce the split_map svg graph [default: F]
        --n             : the file name
        --o             : the output file directory and the median-file directory [default: ./]
        --step          : [all|assembly|others] which steps you want to run Sgrstlfr.
                                all:run the all [default]
                                assembly: only run the module7 to repair the scaffold and genome assembly
                                others: you can appoint which step will run (likes 4,5,6,7) 
        --others        : the others paramaters. likes (minread,10,mdep,200)
        <tools>         : tools path is must !   
USAGE

my ($list,$fasta,$refcover) = @ARGV;
my ($chrtype,$chrid,$psvg,$threads,$aname,$outdir,$step,$others);

GetOptions(
    "chrtype:s" => \$chrtype,
    "chrid:s" => \$chrid,
    "n:s" => \$aname,
    "t:s" => \$threads,
    "others:s" => \$others,
    "p:s" => \$psvg,
    "o:s" => \$outdir,
    "step:s" => \$step
);

die $usage if (@ARGV<3);

$outdir ||=getcwd;
$psvg ||='F';
$aname ||='test';

my $more_option='';
if ($others){
    my @option_set = split/\,/,$others;
    foreach my $pairs (pairs @option_set){
        $more_option .='-'.${$pairs}[0].' '.${$pairs}[1].' ';
    }
}

my $sgrstlfr = "$Bin/run_Sgrstlfr.pl";
my %list;
open INA,$ARGV[0] or die $!;
while (<INA>){
    chomp;
    my ($name,$father,$l1,$l2) = split/\s+/,$_;
    
    $father =~ s/\r\n\s+//g;
    $name =~ s/\r\n\s+//g;
    $l1 =~ s/\r\n\s+//g;
    $l2 =~ s/\r\n\s+//g;
    # p($name);p($father);p($l1);p($l2);exit;
    my $f1 = $father.'/'.$l1;
    my $f2 = $father.'/'.$l2;
    @{$list{$name}} = ($f1,$f2);
}
p(%list);
close INA;
# print Dumper(%list);exit;
open OUTA,">$outdir/$aname.sgrstlfr.sh" or die $!;
print OUTA '#!/bin/bash'."\n";
foreach my $filename (keys %list){
    system("mkdir -p $outdir/$filename");
    my $fastq1 = ${$list{$filename}}[0];
    my $fastq2 = ${$list{$filename}}[1];
    print OUTA "sh $outdir/$filename/$filename.sh \n";
    open OUTB,">$outdir/$filename/$filename.sh" or die $!;
    print OUTB '#!/bin/bash'."\n";
    print OUTB "set -eo pipefail\n";
    print OUTB "echo ==========start at : `date` ==========\n";
    print OUTB "perl $sgrstlfr -fa $fasta -fq1 $fastq1 -fq2 $fastq2 -rfcvg $refcover -n $filename -o $outdir -chrid $chrid -chrtype $chrtype -t $threads -p_svg $psvg -step $step $more_option\n";
    print OUTB "echo ==========end at : `date` ========== && \\\n";
    print OUTB "echo This work is finished!";
    close OUTB;
}
close OUTA;

print "This work is finished!\n";


