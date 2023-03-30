#! /usr/bin/perl -w

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;

Describe: tansforming coverage.....

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 <single cvg> [-option]

      -spid     <str>  : Sample ID/prefix of output file
      -wind     [int]  : length of a window [default: 10]
      -chrid    [str]  : chromosome id
      -outdir          : Directory of output [default: ./]   
      -help            : show this message

USAGE

my ($sample,$wind,$chrid,$outdir,$help);
 
GetOptions(
   "spid:s" => \$sample,
   "wind:s" => \$wind,
   "chrid:s" => \$chrid,
   "outdir:s" => \$outdir,
   "help" => \$help,
);

die $usage if (@ARGV < 1 ||$help);

my $cvg = $ARGV[0];
$outdir ||= getcwd;
$wind ||= 10;

if ($cvg =~ /\.gz$/){
	open IN,"gzip -dc $cvg |" or die $!;
}else{
	open IN,$cvg || die $!;
}

if(defined $ARGV[1]){
  open OUT1,">$ARGV[1]" ||die $!;
}else{
  open OUT1,">$outdir/$sample.tar.depth" ||die $!;
}

if(defined $ARGV[2]){
   open OUT2,">$ARGV[2]" ||die $!;
}else{
   open OUT2,">$outdir/$sample.depth.stat" ||die $!;
}

#=====Finding max value====##
$/= ">";

<IN>;
my @line;

while(<IN>){
	chomp;	
	@line = split /\n/;
	my $chrID = shift @line;
	last if $chrID eq "$chrid";
	}

my @sgdepth;
foreach (@line){
	my @lnsgdep = split /\s+/; 
	push @sgdepth, @lnsgdep;
}
my $numbase = @sgdepth; 

my $maxdep = max (@sgdepth);
print "$numbase\t$maxdep\n";

my %depdist;
my $num = @sgdepth;
my $depth = 0;
my $cycle = int($num/$wind);
foreach (1..$cycle){
   my $n = $_;
   foreach (1..$wind){
     $depth +=	$sgdepth[($n-1)*$wind+$_-1];
   }
   my $depR = $depth/($wind*$maxdep);
   my $averdep = $depth/$wind;
   my $site1 = ($n-1)*$wind+1;
   my $site2 = $site1+$wind-1;
   print OUT1 "$chrid\t$site1\t$site2\t";
   printf OUT1 "%5.4f\t%.2f",$depR,$averdep;
   print OUT1 "\n";
   $depth = 0;
}
my $lstnum = $num -$cycle*$wind;
if($lstnum){

foreach (1..$lstnum){
  $depth +=	$sgdepth[$cycle*$wind+$_-1];
}
my $depR = $depth/($lstnum*$maxdep);
my $averdep = $depth/$lstnum;
my $site = $cycle*$wind+1;  
print OUT1 "$chrid\t$site\t$num\t";
printf OUT1 "%5.4f\t%.2f",$depR,$averdep;
print OUT1 "\n"; 
}

foreach (@sgdepth){
  if(exists $depdist{$_}){ $depdist{$_} ++;}else{$depdist{$_} =1;}  
}

foreach (0..$maxdep){
	 $depdist{$_} =0 unless(exists $depdist{$_});
   my $depdist_percent = ($depdist{$_}/$numbase)*100;
   print OUT2 "$_\t$depdist{$_}\t$depdist_percent\n";
}

close IN;
close OUT1;
close OUT2;

