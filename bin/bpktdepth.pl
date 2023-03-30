#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
#use PerlIO::gzip;
use Cwd;

my $usage=<<USAGE;

Describe: putative recombination event detection by maxdepth
          read type: syn-syn syn-wt wt-wt 
         
  Author: Yun Wang, wangyun\@genomic.cn
 Version: 1.0, Date: 2013-10-13

Usage:perl $0 <rpe.gz/rse.gz>  -[option] 
          
      -n          <str>  : sample ID/prefix of output file [default:test]
      -sp                : sample name
      -o                 : output directory [default: ./]
      -rdstype           : breakpoint supporting read type: rpe or rse [default: rpe] 
      -uniq              : unique mapping allowed
      -t                 : threshold to identify as a breakpoint;
      -sd                : maximal breakpoint coordinate difference allowed
      -help              : show this message
      
USAGE

my ($prefix,$sample,$rdstype,$uniq,$threshold,$sd,$outdir,$help);
GetOptions(
   "n:s" => \$prefix,
   "sp:s" => \$sample,
   "rdstype:s" => \$rdstype,
   "uniq" =>\$uniq,
   "t:s" =>\$threshold,
   "sd:s" =>\$sd,
   "o:s"  => \$outdir,
   "help"  => \$help,
);

die $usage if ($help||@ARGV != 1);

$prefix ||= "test";
$outdir ||= getcwd;
$threshold ||= 3;
$rdstype ||= "rpe";
$sd ||=5 unless ($sd == 0);

if ($ARGV[0] =~ /.gz$/){ #--rpe.gz
  open IN,"gzip -dc $ARGV[0] |" || die $!;
}else{
  open IN, $ARGV[0] ||die $!;
}

open OUTA, ">$outdir/$prefix.bpk.$rdstype.xls" ||die $!;
open OUTB,">$outdir/$prefix.$rdstype.bkpt.log" or die $!;            #JSnum.bkpt.xls

my %recb;
my %stat;
my %rperec;
my $prdsid;

if ($rdstype eq "rpe"){
$stat{"NO"} = 0;
$stat{"REC"} = 0;
$stat{"wt_wt_REC"} = 0;
$stat{"syn_syn_REC"} = 0;
$stat{"wt_syn_REC"} = 0;

while (<IN>){
	my $line1 = $_;
	chomp $line1;
	next unless $line1=~ /\w+/;
	my ($rdsid1,$ori1,$chr1,$mplen1,$site1,$type1,$recom1,$mapnum1) = (split /\s+/,$line1)[0,1,2,3,4,6,8,9];
  $rdsid1 =~ /(.*)\/1$/;
  my $tmpds2 = "$1".'/2';

  my $line2 = <IN>;
  chomp $line2;
  my ($rdsid2,$ori2,$chr2,$mplen2,$site2,$type2,$recom2,$mapnum2) = (split /\s+/,$line2)[0,1,2,3,4,6,8,9];
  print "Warning read info is inconsistent!!! --$rdsid1\t$rdsid2\n"if ($tmpds2 ne $rdsid2 || $type1 ne $type1 || $recom1 ne $recom2);  
  
  my $prdsid_tmp = (split/#/,$rdsid1)[0];
  my $mapnum = $mapnum1 + $mapnum2;
  my $info = "$line1\n$line2\n";
  
  if(defined $uniq){
    next if $mapnum > 2;
  }
  
  if($recom1 eq "NO"){
    $stat{"NO"} ++;
    next;  
  }else{
    $stat{"REC"} ++;
    if($ori1 eq "+"){
       $site1 = $site1 + $mplen1;
    }
    
    if($ori2 eq "-"){
    	 $site2 = $site2 + $mplen2;
    }
    #print "$site1..\n";
    if ($type1 eq "wt_wt" ||$type1 eq "syn_syn"){
 #      
      $stat{"wt_wt_REC"} ++ if $type1 eq "wt_wt";
      $stat{"syn_syn_REC"} ++ if $type1 eq "syn_syn";
      if($chr1 eq $chr2){
        if($site1 < $site2){##
          foreach my $i(-$sd..$sd){
           foreach my $j(-$sd..$sd){
           	   my $sita = $site1 + $i;
           	   my $sitb = $site2 + $j;
           	 
               my $keysid = "$chr1 $sita $chr1 $sitb";
               $recb{$type1}{"DP"}{$keysid}++;
               $recb{$type1}{"INFO"}{$keysid} .= "$info";
               $recb{$type1}{"DIS"}{$keysid} += abs($i)+abs($j);   
           }
         }
        }else{
         foreach my $i(-$sd..$sd){
           foreach my $j(-$sd..$sd){           	
           	   my $sita = $site1 + $i;
           	   my $sitb = $site2 + $j;
               my $keysid = "$chr1 $sitb $chr1 $sita";

               $recb{$type1}{"DP"}{$keysid}++;
               $recb{$type1}{"INFO"}{$keysid} .= "$info";
               $recb{$type1}{"DIS"}{$keysid} += abs($i)+abs($j);     
      	   }
       	}
      
        }
#    
      }else{  ##
       if (ordstr($chr1,1) < ordstr($chr2,1)){
         foreach my $i(-$sd..$sd){
           foreach my $j(-$sd..$sd){
            	my $sita = $site1 + $i;
            	my $sitb = $site2 + $j;
              my $keysid = "$chr1 $sita $chr2 $sitb";

              $recb{$type1}{"DP"}{$keysid}++;
              $recb{$type1}{"INFO"}{$keysid} .= "$info";
              $recb{$type1}{"DIS"}{$keysid} += abs($i)+abs($j);     
         }       	
        }
      }else{
       foreach my $i(-$sd..$sd){
         foreach my $j(-$sd..$sd){
           	my $sita = $site1 + $i;
           	my $sitb = $site2 + $j;
            my $keysid = "$chr2 $sitb $chr1 $sita";
            $recb{$type1}{"DP"}{$keysid}++;
            $recb{$type1}{"INFO"}{$keysid} .= "$info";
            $recb{$type1}{"DIS"}{$keysid} += abs($i)+abs($j);                
         }       	
      }   	
    }
  }      
}elsif($type1 eq "syn_wt"){
     
     $stat{"wt_syn_REC"} ++;
     
     if (ordstr($chr1,1) < ordstr($chr2,1)){
       foreach my $i(-$sd..$sd){
         foreach my $j(-$sd..$sd){
           	my $sita = $site1 + $i;
           	my $sitb = $site2 + $j;
            my $keysid = "$chr1 $sita $chr2 $sitb";

            $recb{$type1}{"DP"}{$keysid}++;
            $recb{$type1}{"INFO"}{$keysid} .= "$info";
            $recb{$type1}{"DIS"}{$keysid} += abs($i)+abs($j);     
         }       	
      }        
   }else{
     foreach my $i(-$sd..$sd){
         foreach my $j(-$sd..$sd){
           	my $sita = $site1 + $i;
           	my $sitb = $site2 + $j;
            my $keysid = "$chr2 $sitb $chr1 $sita";
            $recb{$type1}{"DP"}{$keysid}++;
            $recb{$type1}{"INFO"}{$keysid} .= "$info";        
            $recb{$type1}{"DIS"}{$keysid} += abs($i)+abs($j);        
         }       	
      }   	
    }
  }
 }

}

print OUTB "##$sample\n";
#print OUT "$sample\t$stat{wt_wt_REC}\t$stat{wt_syn_REC}\t$stat{syn_syn_REC}\t";

my %bkpt;
my $tmpbkpt;
#$threshold

foreach my $type("wt_wt","syn_wt","syn_syn"){
  my $maxrdsdepth = 0;
  my $maxrdsdinfo= "";
  
  foreach my $key(keys %{$recb{$type}{"DP"}}){
  	next if $recb{$type}{"DP"}{$key} < $threshold; 
  	 	      
    my ($chr1,$site1,$chr2,$site2) = split /\s+/,$key;  
    my $ck = 0; 
     foreach my $i(-$sd*2..$sd*2){

         foreach my $j(-$sd*2..$sd*2){
           	my $sita = $site1 + $i;
           	my $sitb = $site2 + $j;
            my $keysid = "$chr1 $sita $chr2 $sitb";
            #print ">>$keysid\n";
            if(exists $bkpt{$keysid}){
              $ck ++;
              if ($recb{$type}{"DP"}{$keysid} < $recb{$type}{"DP"}{$key}){
                 delete $bkpt{$keysid};
                 $bkpt{$key} = $recb{$type}{"INFO"}{$key};
              }elsif($recb{$type}{"DP"}{$keysid} == $recb{$type}{"DP"}{$key}&&$recb{$type}{"DIS"}{$keysid}/$recb{$type}{"DP"}{$keysid} > $recb{$type}{"DIS"}{$key}/$recb{$type}{"DP"}{$key}){
                 delete $bkpt{$keysid};
                 $bkpt{$key} = $recb{$type}{"INFO"}{$key};
              }
            }
                     
         }       	
      }
      
     if($ck == 0){
        $bkpt{$key} = $recb{$type}{"INFO"}{$key};
        #print "<<$key\n";    
    }   	        
  }

my %chcbkptid;

###########################
foreach my $key(keys %bkpt){
	    
    my ($chr1,$site1,$chr2,$site2) = split /\s+/,$key;
 
   #CHECK again...
    my $chcrec = 0;
   unless (exists $chcbkptid{$key}){ 
     print OUTA "$chr1\t$site1\t$chr2\t$site2\t$recb{$type}{DP}{$key}\t$type\n";
     print OUTB "#$key\t$recb{$type}{DP}{$key}\n$bkpt{$key}\n";
   }
   
   foreach my $i(-$sd*2..$sd*2){
      foreach my $j(-$sd*2..$sd*2){        
         my $sita = $site1 + $i;
         my $sitb = $site2 + $j;
         my $keysid = "$chr1 $sita $chr2 $sitb";
         $chcbkptid{$keysid} = 1;
       }
    }  
}

undef %bkpt;

}
}elsif($rdstype eq "rse"){
while (<IN>){
	my $line1 = $_;
	chomp $line1;
	my ($rdsid,$ori,$chr,$mplen,$site,$mapnum) = (split /\s+/,$line1)[0,1,2,3,4,9];

  my $prdsid_tmp = (split /#0/,$rdsid)[0];
  my $info = "$line1\n";
  
  if(defined $uniq){
    next if $mapnum > 2;
  }
  
  my $bksite;   
  if(($rdsid =~ /.*\/1$/ && $ori eq "-") || ($rdsid =~ /.*\/2$/ && $ori eq "+")){
    $bksite = $site;
  }elsif(($rdsid =~ /.*\/1$/ && $ori eq "+") || ($rdsid =~ /.*\/2$/ && $ori eq "-")){
    $bksite = $site + $mplen;
  }
  
  foreach my $i(-$sd..$sd){
     my $sita = $bksite + $i;
     my $keysid = "$chr $sita";
     $recb{$keysid}{"DP"}++;
     $recb{$keysid}{"INFO"} .= "$info";
     $recb{$keysid}{"DIS"} += abs($i);   
  }
     
}

print OUTB "##$sample\n";

foreach my $key(keys %recb){
  $recb{$key}{"DIS"} = $recb{$key}{"DIS"}/$recb{$key}{"DP"};
}


my %bkpt;
my $tmpbkpt;
  
foreach my $key(keys %recb){
	
	next if $recb{$key}{"DP"} < $threshold; 
	 	   
  my ($chr,$site);    
  my @arr = split / /,$key;
  $site = pop @arr;
  $chr = join "_",@arr;
  
  my $ck = 0; 
  foreach my $i(-$sd*2..$sd*2){
   	my $sita = $site + $i;
    my $keysid = "$chr\_$sita";
    if(exists $bkpt{$keysid}){
       $ck ++;
       if ($recb{$keysid}{"DP"} < $recb{$key}{"DP"}){
           delete $bkpt{$keysid};
           $bkpt{$key} = $recb{$key}{"INFO"};
       }elsif( $recb{$keysid}{"DP"} == $recb{$key}{"DP"} && $recb{$keysid}{"DIS"} > $recb{$key}{"DIS"}){
           delete $bkpt{$keysid};
           $bkpt{$key} = $recb{$key}{"INFO"}; 
       }	
    }                     
  }       	
      
  if($ck == 0){
     $bkpt{$key} = $recb{$key}{"INFO"};
     #print "<<$key\n";    
 }   	        
}

my %chcbkptid;

#Note split
foreach my $key(keys %bkpt){
    my ($chr,$site);    
    my @arr = split / /,$key;
    $site = pop @arr;
    if(@arr>1){
      $chr = join "_",@arr;
    }else{
    	$chr = $arr[0];
    }
   #CHECK again...

    my $chcrec = 0;
    unless (exists $chcbkptid{$key}){ 
     print OUTA "$chr\t$site\t$recb{$key}{DP}\n";
     print OUTB "#$key\t$recb{$key}{DP}\n$bkpt{$key}\n";
   }
   
   foreach my $i(-$sd*2..$sd*2){     
         my $sita = $site + $i;
         my $keysid = "$chr $sita";
         $chcbkptid{$keysid} = 1;
    }
}  	
	
}

##
sub ordstr{
    my ($chrid,$opero) = @_;
    my $sumasc;
    
    foreach my $str(split //,$chrid){   
      $sumasc += ord($str);
    }
    
    return $sumasc;
}





