#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;

Describe:get the variation reads mapped by bwa...

  Author: Yun Wang, wangyun\@genomic.cn
 Version: 1.0, Date: 2012-11-05

Usage:perl $0 <splitted.map>  [-option]

      --n     <str> : Sample ID/prefix of output file
      --s           : the tools path of samtools   
      --o           : Directory of output [default: ./]   
      --h           : show this message

USAGE

my ($samplename,$outdir,$samtools,$help);
GetOptions(  
   "n:s"   => \$samplename,
   "s:s"  =>\$samtools,
   "o:s" => \$outdir,
   "h"   => \$help
);
die $usage if (@ARGV != 1 ||$help);


$outdir ||= getcwd;
$samtools ||='/hwfssz1/ST_BIOCHEM/P18Z10200N0255/PMO/F14ZQSYJSY1726/pangwending/some_tools/samtools-1.13/bin/samtools';

print "#$samplename\n";
print "#clustering the splitted loxp reads\n";
if ($ARGV[0] =~ /\.sam$/){
  open IN,$ARGV[0] ||die $!;
}else{
  open IN,"$samtools view $ARGV[0] |" or die $!;
}

open OUT1, ">$outdir/$samplename.gvr" ||die $!;
open OUT2, ">$outdir/$samplename.mpi" ||die $!;

my $fline = ();

while(<IN>){
	 chomp;
	 $fline = $_;
	 if($fline =~ /^@/){
	   next;
	 }else{last;}
 }

my ($readid,$ori,$chr,$site,$seq) = (split /\s+/,$fline)[0,1,2,3,9];
my $seqlen = length($seq);
my ($star,$END) =  (0,0);
my $prin = ();
my $prinid = "*";
my $Tnum = 0;

if ($readid =~ /(.+)_\d+$/){$readid = $1;}
if ($ori == 16){                  
	$star  = $site + $seqlen -1;
  $END = $site;
  $ori = "-";                     
 }elsif($ori == 0){
  $star  = $site;
  $END =$site + $seqlen -1;
  $ori = "+";                    
 }elsif($ori == 4){             
  $readid = "*";
 }else{print "There is a unkown flag\n";}

$prin = "$readid\t$chr\t$star\t$END\t$ori";   
while (<IN>){
	chomp;
  my ($readname,$ori2,$chr,$site,$seq) = (split /\s+/,$_)[0,1,2,3,9]; 
  my $seqlen = length($seq);
  my ($star,$END) =  (0,0);
  
 if ($ori2 == 16){
	$star  = $site +$seqlen -1;
  $END = $site;
  $ori2 = "-";
 }elsif($ori2 == 0){
  $star  = $site;
  $END =$site + $seqlen -1;
  $ori2 = "+";
 }elsif($ori2 == 4){
	next;
 }else{print "There is a unkown flag\n";}
 	
 if ($readname =~ /(.+)_\d+$/){$readname = $1;}     
 
 if ($readname eq $readid){
	 
   if ($ori2 eq $ori){print OUT1 "F\t$prin\t$chr\t$star\t$END\t$ori2\n";   
   }else{print OUT1 "R\t$prin\t$chr\t$star\t$END\t$ori2\n";} 
   $Tnum ++;
  }else{
   print OUT2 "$readname\n" if $prinid ne $readname;
   $prinid = $readname;	
  }
   $prin = "$readname\t$chr\t$star\t$END\t$ori2";
	 $readid = $readname;
	 $ori = $ori2;
}
print "#Total number of splited PE reads:$Tnum\n";  

close IN;
close OUT1;
close OUT2;
  
