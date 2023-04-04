#! /usr/bin/perl -w

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
#use PerlIO::gzip;
use Cwd;

my $usage=<<USAGE;

Describe: transfer bpklog to gvr

Author: Yun Wang, wangyun\@genomic.cn
Version: 1.0.0, Date: 2014-11-25

Usage:perl $0 <bkpt.log> -[option]

	-n		: sample name
	-o		: output directory
	-rdsn		: minimal read number to surport a breakpoint [default 2]
	-type		: bpklog type,rpe/rse  [default rpe]
    
USAGE

my $infile = $ARGV[0];

my ($rdsn,$logtype,$name,$outdir);
GetOptions(
	"n:s" => \$name,
	"o:s" => \$outdir,
    "rdsn:s" => \$rdsn,
    "type:s" => \$logtype, #rpe rse
);

$rdsn ||= 2;
$logtype ||= "rpe"; 

die $usage if (!$infile);
open IN,$infile||die $!;
open OUT,">$outdir/$name.all.$logtype.bkpt.gvr" ||die $!;

$/="\n#";
<IN>;
my %prinfo;
while(<IN>){
	chomp;
	my @rdsinfo = split /\n/,$_;
	my $id = shift @rdsinfo;
	
	my @bpkinfo = split /\s+/,$id;
	my $n = pop @bpkinfo;
	next if $n < $rdsn;
	
	if($logtype eq "rpe"){
    my ($chrid,$site);
    
    if ($bpkinfo[0] gt $bpkinfo[2]){
       $chrid = $bpkinfo[2];
       $site =  $bpkinfo[3];
    }elsif($bpkinfo[0] gt $bpkinfo[2]){
       $chrid = $bpkinfo[0];  
        $site = min ($bpkinfo[1],$bpkinfo[3]);
    }else{
       $chrid = $bpkinfo[0];
       $site =  $bpkinfo[1];
    }
    
   foreach my $i(0..$n-1){
	  	my $j = 2*$i;
	  	my @subrdsinfo1 =  split /\s+/,$rdsinfo[$j];
	  	my @subrdsinfo2 =  split /\s+/,$rdsinfo[$j+1];
	    my ($type,$star1,$end1,$star2,$end2);
	  
	    if($subrdsinfo1[1]eq $subrdsinfo2[1]){
	    	 $type = "F";
	    }else{
	    	 $type = "R";
	    }
	  
	    if ($subrdsinfo1[1] eq "+"){
	    	$star1 = $subrdsinfo1[4];
	    	$end1 = $subrdsinfo1[4] + $subrdsinfo1[3];
	    }else{
	   	  $star1 = $subrdsinfo1[4] + $subrdsinfo1[3];
	    	$end1 = $subrdsinfo1[4];
	    }

	    if ($subrdsinfo2[1] eq "+"){
  	  	$star2 = $subrdsinfo2[4];
	    	$end2 = $subrdsinfo2[4] + $subrdsinfo2[3];
	    }else{
	   	  $star2 = $subrdsinfo2[4] + $subrdsinfo2[3];
	    	$end2 = $subrdsinfo2[4];
  	  }
	   if(defined  $chrid){
	     next unless $subrdsinfo1[2]eq $chrid || $subrdsinfo2[2] eq $chrid;
  	 }
	   $prinfo{$chrid}{$site} .= "$type\t$subrdsinfo1[0]\t$subrdsinfo1[2]\t$star1\t$end1\t$subrdsinfo1[1]\t$subrdsinfo2[2]\t$star2\t$end2\t$subrdsinfo2[1]\n";
  	}
  }elsif($logtype eq "rse"){
  	 my ($chrid,$site) = ($bpkinfo[0], $bpkinfo[1]);
  	 
     foreach my $i(0..$n-1){
		  my @subrdsinfo1 =  split /\s+/,$rdsinfo[$i];
		
	    my ($type,$star1,$end1);	 
	  	  $type = "-";
	  
	    if ($subrdsinfo1[1] eq "+"){
	  	  $star1 = $subrdsinfo1[4];
	  	  $end1 = $subrdsinfo1[4] + $subrdsinfo1[3];
	    }else{
	 	    $star1 = $subrdsinfo1[4] + $subrdsinfo1[3];
	  	  $end1 = $subrdsinfo1[4];
	    }
     
     	if (defined $chrid){
		    next if $subrdsinfo1[2] ne $chrid ;
	    }
	   $prinfo{$chrid}{$site} .= "$type\t$subrdsinfo1[0]\t$subrdsinfo1[2]\t$star1\t$end1\t$subrdsinfo1[1]\n";
	  }
  }
}

foreach my $chrid(sort {$a cmp $b} keys %prinfo){
	  foreach my $site(sort {$a <=>$b} keys %{$prinfo{$chrid}}){
	     print OUT "$prinfo{$chrid}{$site}";
	  }

}

close IN;
close OUT;
